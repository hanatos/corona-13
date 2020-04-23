#include "corona_common.h"
#include "colour.h"
#include "display.h"
#include "view.h"
#include "threads.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <xmmintrin.h>
#define XK_LATIN1
#define XK_MISCELLANY
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/keysymdef.h>
#include <X11/XKBlib.h>

typedef struct display_t
{
	int isShuttingDown;
  struct display_t *old_display;
	//Format destFormat;
	Atom wmProtocols;
	Atom wmDeleteWindow;
  Display* display;
  Window window;
  XImage* image;
  GC gc;
  int width;
  int height;
  unsigned int *buffer;   // double buffer to be pushed to x
  float *pixel;           // colour managed output for display
  int bit_depth;
  char msg[256];
  int msg_len;
  int msg_x, msg_y;

  // parallel job state
  uint64_t job_px_counter, job_px_end;
  float job_scale;
  const float *job_buffer;

  // extra controls
  int num_controls;
  int view_controls;
  int active_control;
  display_control_t control[20];

  void (*onKeyDown)(keycode_t);
  void (*onKeyPressed)(keycode_t);
  void (*onKeyUp)(keycode_t);
  void (*onMouseButtonUp)(mouse_t);
  void (*onMouseButtonDown)(mouse_t);
  void (*onMouseMove)(mouse_t);
  void (*onActivate)(char);
  void (*onClose)();
}
display_t;

static const int eventMask = KeyPressMask | KeyReleaseMask | ButtonPressMask | ButtonReleaseMask | PointerMotionMask | ButtonMotionMask;
#define keyMapSize 256

static keycode_t normalKeys[keyMapSize];
static keycode_t functionKeys[keyMapSize];
static char keyIsPressed[keyMapSize];
static char keyIsReleased[keyMapSize];

static int keyMapsInitialized = 0;
static const unsigned char font9x16[] = {
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,30,0,63,176,63,176,30,0,0,0,0,0,0,0,0,0,112,0,120,0,0,0,0,0,120,0,112,0,0,0,0,0,4,64,31,240,31,240,4,64,
4,64,31,240,31,240,4,64,0,0,28,96,62,48,34,16,226,28,226,28,51,240,25,224,0,0,0,0,24,48,24,96,0,192,1,128,3,0,6,0,12,48,24,48,0,0,1,224,27,240,62,16,39,16,61,224,27,240,2,16,0,0,0,0,
0,0,8,0,120,0,112,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,15,192,31,224,48,48,32,16,0,0,0,0,0,0,0,0,0,0,32,16,48,48,31,224,15,192,0,0,0,0,0,0,1,0,5,64,7,192,3,128,3,128,
7,192,5,64,1,0,0,0,0,0,1,0,1,0,7,192,7,192,1,0,1,0,0,0,0,0,0,0,0,0,0,8,0,120,0,112,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,0,0,0,0,0,0,
0,0,0,0,0,48,0,48,0,0,0,0,0,0,0,0,0,48,0,96,0,192,1,128,3,0,6,0,12,0,0,0,0,0,15,192,31,224,48,48,35,16,48,48,31,224,15,192,0,0,0,0,0,0,8,16,24,16,63,240,63,240,0,16,
0,16,0,0,0,0,16,112,48,240,33,144,35,16,38,16,60,48,24,48,0,0,0,0,16,32,48,48,32,16,34,16,34,16,63,240,29,224,0,0,0,0,3,0,7,0,13,0,25,16,63,240,63,240,1,16,0,0,0,0,62,32,62,48,
34,16,34,16,34,16,35,240,33,224,0,0,0,0,15,224,31,240,50,16,34,16,34,16,35,240,1,224,0,0,0,0,48,0,48,0,33,240,35,240,38,0,60,0,56,0,0,0,0,0,29,224,63,240,34,16,34,16,34,16,63,240,29,224,
0,0,0,0,28,0,62,16,34,16,34,16,34,48,63,224,31,192,0,0,0,0,0,0,0,0,0,0,12,96,12,96,0,0,0,0,0,0,0,0,0,0,0,0,0,16,12,112,12,96,0,0,0,0,0,0,0,0,0,0,1,0,3,128,
6,192,12,96,24,48,16,16,0,0,0,0,0,0,4,128,4,128,4,128,4,128,4,128,4,128,0,0,0,0,0,0,16,16,24,48,12,96,6,192,3,128,1,0,0,0,0,0,24,0,56,0,32,0,33,176,35,176,62,0,28,0,0,0,
0,0,15,224,31,240,16,16,19,208,19,208,31,208,15,128,0,0,0,0,7,240,15,240,25,0,49,0,25,0,15,240,7,240,0,0,0,0,32,16,63,240,63,240,34,16,34,16,63,240,29,224,0,0,0,0,15,192,31,224,48,48,32,16,
32,16,48,48,24,96,0,0,0,0,32,16,63,240,63,240,32,16,48,48,31,224,15,192,0,0,0,0,32,16,63,240,63,240,34,16,39,16,48,48,56,112,0,0,0,0,32,16,63,240,63,240,34,16,39,0,48,0,56,0,0,0,0,0,
15,192,31,224,48,48,33,16,33,16,49,224,25,240,0,0,0,0,63,240,63,240,2,0,2,0,2,0,63,240,63,240,0,0,0,0,0,0,0,0,32,16,63,240,63,240,32,16,0,0,0,0,0,0,0,224,0,240,0,16,32,16,63,240,
63,224,32,0,0,0,0,0,32,16,63,240,63,240,3,0,7,128,60,240,56,112,0,0,0,0,32,16,63,240,63,240,32,16,0,16,0,48,0,112,0,0,0,0,63,240,63,240,28,0,14,0,28,0,63,240,63,240,0,0,0,0,63,240,
63,240,28,0,14,0,7,0,63,240,63,240,0,0,0,0,31,224,63,240,32,16,32,16,32,16,63,240,31,224,0,0,0,0,32,16,63,240,63,240,34,16,34,0,62,0,28,0,0,0,0,0,31,224,63,240,32,16,32,48,32,28,63,252,
31,228,0,0,0,0,32,16,63,240,63,240,34,0,35,0,63,240,28,240,0,0,0,0,24,96,60,112,38,16,34,16,35,16,57,240,24,224,0,0,0,0,0,0,56,0,48,16,63,240,63,240,48,16,56,0,0,0,0,0,63,224,63,240,
0,16,0,16,0,16,63,240,63,224,0,0,0,0,63,128,63,192,0,96,0,48,0,96,63,192,63,128,0,0,0,0,63,224,63,240,0,112,3,192,0,112,63,240,63,224,0,0,0,0,48,48,60,240,15,192,7,128,15,192,60,240,48,48,
0,0,0,0,0,0,60,0,62,16,3,240,3,240,62,16,60,0,0,0,0,0,56,112,48,240,33,144,35,16,38,16,60,48,56,112,0,0,0,0,0,0,0,0,63,240,63,240,32,16,32,16,0,0,0,0,0,0,28,0,14,0,7,0,
3,128,1,192,0,224,0,112,0,0,0,0,0,0,0,0,32,16,32,16,63,240,63,240,0,0,0,0,0,0,16,0,48,0,96,0,192,0,96,0,48,0,16,0,0,0,0,0,0,4,0,4,0,4,0,4,0,4,0,4,0,4,0,4,
0,0,0,0,0,0,96,0,112,0,16,0,0,0,0,0,0,0,0,0,0,224,5,240,5,16,5,16,7,224,3,240,0,16,0,0,0,0,32,16,63,240,63,224,4,16,6,16,3,240,1,224,0,0,0,0,3,224,7,240,4,16,4,16,
4,16,6,48,2,32,0,0,0,0,1,224,3,240,6,16,36,16,63,224,63,240,0,16,0,0,0,0,3,224,7,240,5,16,5,16,5,16,7,48,3,32,0,0,0,0,0,0,2,16,31,240,63,240,34,16,48,0,24,0,0,0,0,0,
3,228,7,246,4,18,4,18,3,254,7,252,4,0,0,0,0,0,32,16,63,240,63,240,2,0,4,0,7,240,3,240,0,0,0,0,0,0,0,0,4,16,55,240,55,240,0,16,0,0,0,0,0,0,0,0,0,4,0,6,0,2,4,2,
55,254,55,252,0,0,0,0,32,16,63,240,63,240,1,128,3,192,6,112,4,48,0,0,0,0,0,0,0,0,32,16,63,240,63,240,0,16,0,0,0,0,0,0,7,240,7,240,6,0,3,240,3,240,6,0,7,240,3,240,0,0,4,0,
7,240,3,240,4,0,4,0,7,240,3,240,0,0,0,0,3,224,7,240,4,16,4,16,4,16,7,240,3,224,0,0,0,0,4,2,7,254,3,254,4,18,4,16,7,240,3,224,0,0,0,0,3,224,7,240,4,16,4,18,3,254,7,254,
4,2,0,0,0,0,4,16,7,240,3,240,6,16,4,0,6,0,2,0,0,0,0,0,3,32,7,176,4,144,4,144,4,144,6,240,2,96,0,0,0,0,4,0,4,0,31,224,63,240,4,16,4,48,0,32,0,0,0,0,7,224,7,240,
0,16,0,16,7,224,7,240,0,16,0,0,0,0,7,192,7,224,0,48,0,16,0,48,7,224,7,192,0,0,0,0,7,224,7,240,0,48,0,224,0,224,0,48,7,240,7,224,0,0,4,16,6,48,3,96,1,192,1,192,3,96,6,48,
4,16,0,0,7,226,7,242,0,18,0,18,0,22,7,252,7,248,0,0,0,0,6,48,6,112,4,208,5,144,7,16,6,48,4,48,0,0,0,0,0,0,2,0,2,0,31,224,61,240,32,16,32,16,0,0,0,0,0,0,0,0,0,0,
62,248,62,248,0,0,0,0,0,0,0,0,0,0,32,16,32,16,61,240,31,224,2,0,2,0,0,0,0,0,32,0,96,0,64,0,96,0,32,0,96,0,64,0,0,0,0,0,1,224,3,224,6,32,12,32,6,32,3,224,1,224,0,0};

int initializeKeyMaps()
{
  for (int i = 0; i < keyMapSize; ++i)
  {
    normalKeys[i] = KeyUndefined;
    functionKeys[i] = KeyUndefined;
    keyIsPressed[i] = 0;
    keyIsReleased[i] = 0;
  }

  normalKeys[XK_space] = KeySpace;
  normalKeys[XK_comma] = KeyComma;
  normalKeys[XK_period] = KeyPeriod;
  normalKeys[XK_slash] = KeySlash;
  normalKeys[XK_0] = KeyZero;
  normalKeys[XK_1] = KeyOne;
  normalKeys[XK_2] = KeyTwo;
  normalKeys[XK_3] = KeyThree;
  normalKeys[XK_4] = KeyFour;
  normalKeys[XK_5] = KeyFive;
  normalKeys[XK_6] = KeySix;
  normalKeys[XK_7] = KeySeven;
  normalKeys[XK_8] = KeyEight;
  normalKeys[XK_9] = KeyNine;
  normalKeys[XK_semicolon] = KeySemiColon;
  normalKeys[XK_equal] = KeyEquals;
  normalKeys[XK_a] = KeyA;
  normalKeys[XK_b] = KeyB;
  normalKeys[XK_c] = KeyC;
  normalKeys[XK_d] = KeyD;
  normalKeys[XK_e] = KeyE;
  normalKeys[XK_f] = KeyF;
  normalKeys[XK_g] = KeyG;
  normalKeys[XK_h] = KeyH;
  normalKeys[XK_i] = KeyI;
  normalKeys[XK_j] = KeyJ;
  normalKeys[XK_k] = KeyK;
  normalKeys[XK_l] = KeyL;
  normalKeys[XK_m] = KeyM;
  normalKeys[XK_n] = KeyN;
  normalKeys[XK_o] = KeyO;
  normalKeys[XK_p] = KeyP;
  normalKeys[XK_q] = KeyQ;
  normalKeys[XK_r] = KeyR;
  normalKeys[XK_s] = KeyS;
  normalKeys[XK_t] = KeyT;
  normalKeys[XK_u] = KeyU;
  normalKeys[XK_v] = KeyV;
  normalKeys[XK_w] = KeyW;
  normalKeys[XK_x] = KeyX;
  normalKeys[XK_y] = KeyY;
  normalKeys[XK_z] = KeyZ;
  normalKeys[XK_bracketleft] = KeyOpenBracket;
  normalKeys[XK_backslash] = KeyBackSlash;
  normalKeys[XK_bracketright] = KeyCloseBracket;

  functionKeys[0xff & XK_BackSpace] = KeyBackSpace;
  functionKeys[0xff & XK_Tab] = KeyTab;
  functionKeys[0xff & XK_Linefeed] = KeyUndefined;
  functionKeys[0xff & XK_Clear] = KeyClear;
  functionKeys[0xff & XK_Return] = KeyEnter;
  functionKeys[0xff & XK_Pause] = KeyPause;
  functionKeys[0xff & XK_Scroll_Lock] = KeyScrollLock;
  functionKeys[0xff & XK_Sys_Req] = KeyPrintScreen;
  functionKeys[0xff & XK_Escape] = KeyEscape;
  functionKeys[0xff & XK_Delete] = KeyDelete;
  functionKeys[0xff & XK_Kanji] = KeyKanji;
  functionKeys[0xff & XK_Kana_Shift] = KeyKana;
  functionKeys[0xff & XK_Home] = KeyHome;
  functionKeys[0xff & XK_Left] = KeyLeft;
  functionKeys[0xff & XK_Up] = KeyUp;
  functionKeys[0xff & XK_Right] = KeyRight;
  functionKeys[0xff & XK_Down] = KeyDown;
  functionKeys[0xff & XK_Prior] = KeyUndefined;
  functionKeys[0xff & XK_Page_Up] = KeyPageUp;
  functionKeys[0xff & XK_Next] = KeyUndefined;
  functionKeys[0xff & XK_Page_Down] = KeyPageDown;
  functionKeys[0xff & XK_End] = KeyEnd;
  functionKeys[0xff & XK_Begin] = KeyUndefined;
  functionKeys[0xff & XK_Select] = KeyUndefined;
  functionKeys[0xff & XK_Print] = KeyUndefined;
  functionKeys[0xff & XK_Execute] = KeyUndefined;
  functionKeys[0xff & XK_Insert] = KeyInsert;
  functionKeys[0xff & XK_Undo] = KeyUndefined;
  functionKeys[0xff & XK_Redo] = KeyUndefined;
  functionKeys[0xff & XK_Menu] = KeyUndefined;
  functionKeys[0xff & XK_Find] = KeyUndefined;
  functionKeys[0xff & XK_Cancel] = KeyCancel;
  functionKeys[0xff & XK_Help] = KeyHelp;
  functionKeys[0xff & XK_Break] = KeyUndefined;
  functionKeys[0xff & XK_Mode_switch] = KeyModeChange;
  functionKeys[0xff & XK_Num_Lock] = KeyNumLock;
  functionKeys[0xff & XK_KP_Space] = KeySpace;
  functionKeys[0xff & XK_KP_Tab] = KeyTab;
  functionKeys[0xff & XK_KP_Enter] = KeyEnter;
  functionKeys[0xff & XK_KP_F1] = KeyF1;
  functionKeys[0xff & XK_KP_F2] = KeyF2;
  functionKeys[0xff & XK_KP_F3] = KeyF3;
  functionKeys[0xff & XK_KP_F4] = KeyF4;
  functionKeys[0xff & XK_KP_Home] = KeyHome;
  functionKeys[0xff & XK_KP_Left] = KeyLeft;
  functionKeys[0xff & XK_KP_Right] = KeyRight;
  functionKeys[0xff & XK_KP_Down] = KeyDown;
  functionKeys[0xff & XK_KP_Prior] = KeyUndefined;
  functionKeys[0xff & XK_KP_Page_Up] = KeyPageUp;
  functionKeys[0xff & XK_KP_Next] = KeyUndefined;
  functionKeys[0xff & XK_KP_Page_Down] = KeyPageDown;
  functionKeys[0xff & XK_KP_End] = KeyEnd;
  functionKeys[0xff & XK_KP_Begin] = KeyUndefined;
  functionKeys[0xff & XK_KP_Insert] = KeyInsert;
  functionKeys[0xff & XK_KP_Delete] = KeyDelete;
  functionKeys[0xff & XK_KP_Equal] = KeyEquals;
  functionKeys[0xff & XK_KP_Multiply] = KeyMultiply;
  functionKeys[0xff & XK_KP_Add] = KeyAdd;
  functionKeys[0xff & XK_KP_Separator] = KeySeparator;
  functionKeys[0xff & XK_KP_Subtract] = KeySubtract;
  functionKeys[0xff & XK_KP_Decimal] = KeyDecimal;
  functionKeys[0xff & XK_KP_Divide] = KeyDivide;
  functionKeys[0xff & XK_KP_0] = KeyNumPad0;
  functionKeys[0xff & XK_KP_1] = KeyNumPad1;
  functionKeys[0xff & XK_KP_2] = KeyNumPad2;
  functionKeys[0xff & XK_KP_3] = KeyNumPad3;
  functionKeys[0xff & XK_KP_4] = KeyNumPad4;
  functionKeys[0xff & XK_KP_5] = KeyNumPad5;
  functionKeys[0xff & XK_KP_6] = KeyNumPad6;
  functionKeys[0xff & XK_KP_7] = KeyNumPad7;
  functionKeys[0xff & XK_KP_8] = KeyNumPad8;
  functionKeys[0xff & XK_KP_9] = KeyNumPad9;
  functionKeys[0xff & XK_F1] = KeyF1;
  functionKeys[0xff & XK_F2] = KeyF2;
  functionKeys[0xff & XK_F3] = KeyF3;
  functionKeys[0xff & XK_F4] = KeyF4;
  functionKeys[0xff & XK_F5] = KeyF5;
  functionKeys[0xff & XK_F6] = KeyF6;
  functionKeys[0xff & XK_F7] = KeyF7;
  functionKeys[0xff & XK_F8] = KeyF8;
  functionKeys[0xff & XK_F9] = KeyF9;
  functionKeys[0xff & XK_F10] = KeyF10;
  functionKeys[0xff & XK_F11] = KeyF11;
  functionKeys[0xff & XK_F12] = KeyF12;
  functionKeys[0xff & XK_Shift_L] = KeyShift;
  functionKeys[0xff & XK_Shift_R] = KeyShift;
  functionKeys[0xff & XK_Control_L] = KeyControl;
  functionKeys[0xff & XK_Control_R] = KeyControl;
  functionKeys[0xff & XK_Caps_Lock] = KeyCapsLock;
  functionKeys[0xff & XK_Shift_Lock] = KeyCapsLock;
  functionKeys[0xff & XK_Meta_L] = KeyMeta;
  functionKeys[0xff & XK_Meta_R] = KeyMeta;
  functionKeys[0xff & XK_Alt_L] = KeyAlt;
  functionKeys[0xff & XK_Alt_R] = KeyAlt;
  return 1;
}

display_t *display_open(const char title[], int width, int height)
{
  if(!keyMapsInitialized) keyMapsInitialized = initializeKeyMaps();
  display_t *d = (display_t*) malloc(sizeof(display_t));

  d->num_controls = 0;
  d->view_controls = 0;
  d->active_control = 0;

  d->old_display = NULL;
  // let's open a display

  d->display = XOpenDisplay(0);
  if (!d->display)
    return 0;

  d->width = width;
  d->height = height;
  d->isShuttingDown = 0;
  d->onKeyDown = NULL;
  d->onKeyPressed = NULL;
  d->onKeyUp = NULL;
  d->onMouseButtonUp = NULL;
  d->onMouseButtonDown = NULL;
  d->onMouseMove = NULL;
  d->onActivate = NULL;
  d->onClose = NULL;

  const int screen = DefaultScreen(d->display);
  Visual* visual = DefaultVisual(d->display, screen);
  if (!visual)
  {
    display_close(d);
    return 0;
  }

  // It gets messy when talking about color depths.
  //
  // For the image buffer, we either need 8, 16 or 32 bitsPerPixel.  8 bits we'll 
  // never have (hopefully), 16 bits will be used for displayDepth 15 & 16, and 
  // 32 bits must be used for depths 24 and 32.
  //
  // The converters will get this right when talking about displayDepth 15 & 16, but 
  // it will wrongly assume that displayDepth 24 takes _exactly_ 24 bitsPerPixel.  We 
  // solve that by tricking the converter requester by presenting it a 32 bit
  // bufferDepth instead.
  //
  const int displayDepth = DefaultDepth(d->display, screen);
  d->bit_depth = displayDepth;
  const int bufferDepth = displayDepth == 24 ? 32 : displayDepth;
  const int bytesPerPixel = (bufferDepth + 7) / 8;
  const int bitsPerPixel = 8 * bytesPerPixel;
  if (bitsPerPixel != 16 && bitsPerPixel != 32)
  {
    display_close(d);
    return 0;
  }
#if 0
		destFormat_ = findFormat(bufferDepth,
			visual->red_mask, visual->green_mask, visual->blue_mask);
		floatingPointConverter_ = requestConverter(Format::XBGRFFFF, destFormat_);
		trueColorConverter_ = requestConverter(Format::XRGB8888, destFormat_);
		if (!floatingPointConverter_ || !trueColorConverter_)
		{
			close();
			return false;
		}
#endif

  // let's create a window
		
  const Window root = DefaultRootWindow(d->display);
		
  const int screenWidth = DisplayWidth(d->display, screen);
  const int screenHeight = DisplayHeight(d->display, screen);
  const int left = (screenWidth - width) / 2;
  const int top = (screenHeight - height) / 2;

  XSetWindowAttributes attributes;
  attributes.border_pixel = attributes.background_pixel = BlackPixel(d->display, screen);
  attributes.backing_store = NotUseful;

  d->window = XCreateWindow(d->display, root, left, top, width, height, 0,
      displayDepth, InputOutput, visual, 
      CWBackPixel | CWBorderPixel | CWBackingStore, &attributes);


  XStoreName(d->display, d->window, title);

  d->wmProtocols = XInternAtom(d->display, "WM_PROTOCOLS", True);
  d->wmDeleteWindow = XInternAtom(d->display, "WM_DELETE_WINDOW", True);
  if (d->wmProtocols == 0 || d->wmDeleteWindow == 0)
  {
    display_close(d);
    return 0;
  }
  if (XSetWMProtocols(d->display, d->window, &(d->wmDeleteWindow), 1) == 0)
  {
    display_close(d);
    return 0;
  }

  XSizeHints sizeHints;
  sizeHints.flags = PPosition | PMinSize | PMaxSize;
  sizeHints.x = sizeHints.y = 0;
  sizeHints.min_width = sizeHints.max_width = width;
  sizeHints.min_height = sizeHints.max_height = height;
  XSetNormalHints(d->display, d->window, &sizeHints);
  XClearWindow(d->display, d->window);
  XSelectInput(d->display, d->window, eventMask);

  // create (image) buffer

  d->pixel = common_alloc(16, sizeof(float)*width*height*4);
  memset(d->pixel, 0, sizeof(float)*width*height*4);

  d->buffer = (unsigned int*) common_alloc(128, sizeof(char) * width * height * bytesPerPixel);
  if (!d->buffer)
  {
    display_close(d);
    return 0;
  }

  d->gc = DefaultGC(d->display, screen);
  d->image = XCreateImage(d->display, CopyFromParent, displayDepth, ZPixmap, 0, 0,
      width, height, bitsPerPixel, width * bytesPerPixel);
#if 1//defined(__LITTLE_ENDIAN__)
  d->image->byte_order = LSBFirst;
#else
  d->image->byte_order = MSBFirst;
#endif	
  if (!d->image)
  {
    display_close(d);
    return 0;
  }

  d->msg[0] = '\0';
  d->msg_len = 0;

  // we have a winner!

  XMapRaised(d->display, d->window);
  XFlush(d->display);

  return d;
}

void display_close(display_t *d)
{	
  if (d->image)
    XDestroyImage(d->image);

  if (d->display && d->window)
    XDestroyWindow(d->display, d->window);

  if (d->display)
    XCloseDisplay(d->display);

  free(d->buffer);
  free(d->pixel);
  free(d);
}

static inline void display_control_update(
    display_t *d,
    int num,
    int inc)
{
  if(num < 0 || num >= d->num_controls) return;
  if(d->control[num].logscale)
  {
    if(inc) d->control[num].storage[0] = CLAMP(d->control[num].storage[0] * d->control[num].step, d->control[num].min, d->control[num].max);
    else    d->control[num].storage[0] = CLAMP(d->control[num].storage[0] / d->control[num].step, d->control[num].min, d->control[num].max);
  }
  else
  {
    if(inc) d->control[num].storage[0] = CLAMP(d->control[num].storage[0] + d->control[num].step, d->control[num].min, d->control[num].max);
    else    d->control[num].storage[0] = CLAMP(d->control[num].storage[0] - d->control[num].step, d->control[num].min, d->control[num].max);
    // round to nearest step
    d->control[num].storage[0] = ((int)(d->control[num].storage[0] / d->control[num].step + .5f)) * d->control[num].step;
  }
  if(d->control[num].clear) view_clear();
}

static inline void handleEvent(const XEvent *event, display_t *d)
{
  switch (event->type)
  {
    case KeyPress:
    case KeyRelease:
      {
        const KeySym keySym = XkbKeycodeToKeysym (d->display, event->xkey.keycode, 0, 0);
        const int hiSym = (keySym & 0xff00) >> 8;
        const int loSym = keySym & 0xff;

        keycode_t code = KeyUndefined;
        switch (hiSym)
        {
          case 0x00:
            code = normalKeys[loSym];
            break;
          case 0xff:
            code = functionKeys[loSym];
            break;
        }

        if (event->type == KeyPress)
        {
          if (code == KeyTab)
          {
            d->view_controls = 1-d->view_controls;
            break;
          }
          if(d->view_controls && code == KeyDown && d->active_control < d->num_controls-1)
            d->active_control++;
          else if(d->view_controls && code == KeyUp && d->active_control > 0)
            d->active_control--;
          else if(d->view_controls && code == KeyLeft)
            display_control_update(d, d->active_control, 0);
          else if(d->view_controls && code == KeyRight)
            display_control_update(d, d->active_control, 1);
          else if (!keyIsPressed[code])
          {
            if(d->onKeyDown) d->onKeyDown(code);
            else if (code == KeyEscape) d->isShuttingDown = 1;
          }
          keyIsPressed[code] = 1;
          keyIsReleased[code] = 0;
        }
        else
        {
          keyIsReleased[code] = 1;
        }
        break;
      }

    case ButtonPress:
    case ButtonRelease:
      {
        mouse_t mouse;
        mouse.x = event->xbutton.x;
        mouse.y = event->xbutton.y;
        mouse.buttons.left = event->xbutton.button == Button1;
        mouse.buttons.middle = event->xbutton.button == Button2;
        mouse.buttons.right = event->xbutton.button == Button3;
        if (event->type == ButtonPress)
        {
          if (d->onMouseButtonDown) d->onMouseButtonDown(mouse);
        }
        else
        {
          if (d->onMouseButtonUp) d->onMouseButtonUp(mouse);
        }
        break;
      }
    case MotionNotify:
      {
        mouse_t mouse;
        mouse.x = event->xmotion.x;
        mouse.y = event->xmotion.y;
        mouse.buttons.left = (event->xmotion.state & Button1Mask) != 0;
        mouse.buttons.middle = (event->xmotion.state & Button2Mask) != 0;
        mouse.buttons.right = (event->xmotion.state & Button3Mask) != 0;
        if (d->onMouseMove) d->onMouseMove(mouse);
        break;
      }
    case ClientMessage:
      {
        if (event->xclient.message_type == d->wmProtocols && 
            event->xclient.format == 32 &&
            event->xclient.data.l[0] == (long) d->wmDeleteWindow)
        {
          if(d->onClose) d->onClose();
          else d->isShuttingDown = 1;
        }
        break;
      }
  }
}

void display_pump_events(display_t *d)
{
  /*if(d->old_display)
  {
    display_close(d->old_display);
    d->old_display = NULL;
  }*/
  XEvent event;
  while (1)
  {		
    //if (d->isShuttingDown) return;
    if (XCheckWindowEvent(d->display, d->window, -1, &event))
      handleEvent(&event, d);
    else if (XCheckTypedEvent(d->display, ClientMessage, &event))
      handleEvent(&event, d);
    else break;
  }

  // send key press and up events

  for (int i = 0; i < keyMapSize; ++i)
  {
    if (keyIsReleased[i] && keyIsPressed[i])
    {
      if (d->onKeyUp) d->onKeyUp((keycode_t)i);
      keyIsPressed[i] = 0;
      keyIsReleased[i] = 0;
    }
    else if (keyIsPressed[i])
    {
      if (d->onKeyPressed) d->onKeyPressed((keycode_t)i);
    }
  }		
}

static inline void display_render_text(
    display_t *d,
    const char *text,
    int text_x,
    int text_y)
{
  int px = text_x;
  for (int pos = 0; text[pos] != '\0'; pos++)
  {
    int charPos = (text[pos] - 32)*9*2;
    for (int x = 0; x < 9*2; x++)
    {
      if (px >= d->width) return;
      unsigned char cLine = font9x16[charPos+x];			
      for (int i = 0; i < 8; i++)
      {
        int y = text_y - (15 - (i + (x&1)*8));
        if ((y >= d->height) || (y < 0)) return;
        if (cLine & (1<<(7-i)))
          d->buffer[px + y*d->width] = -1;
        else
          d->buffer[px + y*d->width] = 0;
      }
      if (x&1) px++;
    }
  }
}

static void *colourmanage(void *arg)
{
  display_t *d = arg;
  const float inv_o = d->job_scale;
  while(1)
  {
    uint64_t i = __sync_fetch_and_add(&d->job_px_counter, 1);
    if(i >= d->job_px_end) return 0;
    // convert camera/framebuffer to display profile:
    float xyz[3], rgb[3];
    for(int k=0;k<3;k++) rgb[k] = d->job_buffer[3*i+k]*inv_o;

    colour_camera_to_xyz(rgb, xyz);
    colour_xyz_to_output(xyz, d->pixel + 4*i + 1);
  }
}

int display_update(display_t *d, const float *input, const float scale)
{
  if (d->isShuttingDown)
  {
    display_close(d);
    return 0;
  }

  if (!d->display || !d->window || !d->image)
    return 0;

  // colour manage for display
  threads_t *t = rt.threads;
  // loop over all pixels once
  d->job_px_counter = 0;
  d->job_px_end = view_width()*view_height();
  d->job_scale = scale;
  d->job_buffer = input;
  for(int k=0;k<rt.num_threads;k++)
    pthread_pool_task_init(t->task + k, &t->pool, colourmanage, d);
  pthread_pool_wait(&t->pool);

  const int w = d->width;
  const int h = d->height;
  const int size = w * h;
  int index = 0;
  const float *pixels = d->pixel;
  if(d->bit_depth == 30)
  {
    for (int i = 0; i < size; i++)
    {
      // test 10-bit vs 8-bit gradient:
      // const int r = (int)(0x3ff*CLAMP((i%w)/(w-1.0), 0, 1)) << 20;
      // const int g = (int)(0x3ff*CLAMP((i%w)/(w-1.0), 0, 1)) << 10;
      // const int b = (int)(0x3ff*CLAMP((i%w)/(w-1.0), 0, 1)) << 0;
      const int r = (int)(0x3ff*CLAMP(pixels[index+3], 0, 1)) << 20;
      const int g = (int)(0x3ff*CLAMP(pixels[index+2], 0, 1)) << 10;
      const int b = (int)(0x3ff*CLAMP(pixels[index+1], 0, 1));
      index += 4;
      // rgba
      d->buffer[i] = r | g | b;
    }
  }
  else // assume 24
  {
    for (int i = 0; i < size; i++)
    {
      const int r = (int)(0xff*CLAMP(pixels[index+1], 0, 1)) << 16;
      const int g = (int)(0xff*CLAMP(pixels[index+2], 0, 1)) << 8;
      const int b = (int)(0xff*CLAMP(pixels[index+3], 0, 1));
      index += 4;
      // rgba
      d->buffer[i] = r | g | b;
    }
  }
  // render message:
  display_render_text(d, d->msg, d->msg_x, d->msg_y);

  // render controls:
  if(d->view_controls) for(int c=0;c<d->num_controls;c++)
  {
    char text[40];
    snprintf(text, 40, "%c %s %g", c==d->active_control ? '*' : ' ', d->control[c].name, *d->control[c].storage);
    display_render_text(d, text, d->width-25*8, 15 + 16*c);
  }

  d->image->data = (char*)(d->buffer);

  XPutImage(d->display, d->window, d->gc, d->image, 0, 0, 0, 0, w, h);
  XFlush(d->display);

  d->image->data = NULL;
  return 1;
}

int display_update_rgba(display_t *d, const unsigned int * rgba)
{
  if (d->isShuttingDown)
  {
    display_close(d);
    return 0;
  }

  if (!d->display || !d->window || !d->image)
    return 0;

  const int w = d->width;
  const int h = d->height;

  if(d->bit_depth == 30)
  {
    for (int i = 0; i < w*h; i++)
    {
      const int r = (((rgba[i]>>16)&0xff)<<2);
      const int g = (((rgba[i]>> 8)&0xff)<<2) << 10;
      const int b = (((rgba[i]    )&0xff)<<2) << 20;
      // 30-bit rgba
      d->buffer[i] = r | g | b;
    }
    d->image->data = (char*)(d->buffer);
  }
  else
  { // assume 24
    d->image->data = (char*)rgba;
  }

  XPutImage(d->display, d->window, d->gc, d->image, 0, 0, 0, 0, w, h);
  XFlush(d->display);

  d->image->data = NULL;
  return 1;
}

display_t *display_resize(display_t *d, const char *title, const int w, const int h)
{
#if 0
  //XResizeWindow(d->display, d->window, w, h);
  display_t *d3 = display_open(title, w, h);
  d3->onKeyDown = d->onKeyDown;
  d3->onKeyUp = d->onKeyUp;
  d3->onMouseButtonDown = d->onMouseButtonDown;
  d3->onMouseMove = d->onMouseMove;
  d3->onClose = d->onClose;
  d->isShuttingDown = 1;
  //display_close(d); // <- segfault.
  d3->old_display = d;
#endif
  return d;
}

void display_print(display_t *d, const int px, const int py, const char *msg, ...)
{
  va_list ap;
  va_start(ap, msg);
  vsnprintf(d->msg, 255, msg, ap);
  // vprintf(msg, ap);
  // printf("\n");
  va_end(ap);
  d->msg_len = strlen(d->msg);
  d->msg_x = px;
  d->msg_y = py + 15;
}

int display_control_add(display_t *d, const char *name, float *storage, float min, float max, float step, int logscale, int clear)
{
  if(d->num_controls >= 20) return 1;
  display_control_t *c = d->control + d->num_controls;
  strncpy(c->name, name, 20);
  c->storage = storage;
  c->min = min;
  c->max = max;
  c->step = step;
  c->clear = clear;
  c->logscale = logscale;
  d->num_controls++;
  return 0;
}

void display_print_info(FILE *fd)
{
  fprintf(fd, "display  : xorg\n");
}

void display_print_usage() { }

void display_register_callbacks(display_t *d, 
  void (*onKeyDown)(keycode_t),
  void (*onKeyPressed)(keycode_t),
  void (*onKeyUp)(keycode_t),
  void (*onMouseButtonDown)(mouse_t),
  void (*onMouseButtonUp)(mouse_t),
  void (*onMouseMove)(mouse_t),
  void (*onActivate)(char),
  void (*onClose)())
{
  d->onKeyDown = onKeyDown;
  d->onKeyPressed = onKeyPressed;
  d->onKeyUp = onKeyUp;
  d->onMouseButtonDown = onMouseButtonDown;
  d->onMouseButtonUp = onMouseButtonUp;
  d->onMouseMove = onMouseMove;
  d->onActivate = onActivate;
  d->onClose = onClose;
}
