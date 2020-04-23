/*
    This file is part of corona-13.

    corona-13 is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    corona-13 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with corona-13. If not, see <http://www.gnu.org/licenses/>.
*/

#include "display.h"
#include "corona_common.h"
#include <stdlib.h>
#include <GL/glew.h>
#include <GL/gl.h>
#include <SDL/SDL.h>



#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#ifdef __SSE__
  #include <xmmintrin.h>
#endif

#define keyMapSize 256

keycode_t normalKeys[keyMapSize];
// keycode_t functionKeys[keyMapSize];
// char keyIsPressed[keyMapSize];
// char keyIsReleased[keyMapSize];

int keyMapsInitialized = 0;
const unsigned char font9x16[] = {
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
    // functionKeys[i] = KeyUndefined;
  }

  normalKeys[SDLK_SPACE] = KeySpace;
  normalKeys[SDLK_COMMA] = KeyComma;
  normalKeys[SDLK_PERIOD] = KeyPeriod;
  normalKeys[SDLK_SLASH] = KeySlash;
  normalKeys[SDLK_0] = KeyZero;
  normalKeys[SDLK_1] = KeyOne;
  normalKeys[SDLK_2] = KeyTwo;
  normalKeys[SDLK_3] = KeyThree;
  normalKeys[SDLK_4] = KeyFour;
  normalKeys[SDLK_5] = KeyFive;
  normalKeys[SDLK_6] = KeySix;
  normalKeys[SDLK_7] = KeySeven;
  normalKeys[SDLK_8] = KeyEight;
  normalKeys[SDLK_9] = KeyNine;
  normalKeys[SDLK_SEMICOLON] = KeySemiColon;
  normalKeys[SDLK_EQUALS] = KeyEquals;
  normalKeys[SDLK_a] = KeyA;
  normalKeys[SDLK_b] = KeyB;
  normalKeys[SDLK_c] = KeyC;
  normalKeys[SDLK_d] = KeyD;
  normalKeys[SDLK_e] = KeyE;
  normalKeys[SDLK_f] = KeyF;
  normalKeys[SDLK_g] = KeyG;
  normalKeys[SDLK_h] = KeyH;
  normalKeys[SDLK_i] = KeyI;
  normalKeys[SDLK_j] = KeyJ;
  normalKeys[SDLK_k] = KeyK;
  normalKeys[SDLK_l] = KeyL;
  normalKeys[SDLK_m] = KeyM;
  normalKeys[SDLK_n] = KeyN;
  normalKeys[SDLK_o] = KeyO;
  normalKeys[SDLK_p] = KeyP;
  normalKeys[SDLK_q] = KeyQ;
  normalKeys[SDLK_r] = KeyR;
  normalKeys[SDLK_s] = KeyS;
  normalKeys[SDLK_t] = KeyT;
  normalKeys[SDLK_u] = KeyU;
  normalKeys[SDLK_v] = KeyV;
  normalKeys[SDLK_w] = KeyW;
  normalKeys[SDLK_x] = KeyX;
  normalKeys[SDLK_y] = KeyY;
  normalKeys[SDLK_z] = KeyZ;
  normalKeys[SDLK_LEFTBRACKET] = KeyOpenBracket;
  normalKeys[SDLK_BACKSLASH] = KeyBackSlash;
  normalKeys[SDLK_RIGHTBRACKET] = KeyCloseBracket;
  return 1;
}

display_t *display_open(const char title[], int width, int height)
{
  if(!keyMapsInitialized) keyMapsInitialized = initializeKeyMaps();
  display_t *d = (display_t*) malloc(sizeof(display_t));

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
  d->msg[0] = '\0';
  d->msg_len = 0;

  const SDL_VideoInfo* info = NULL;
  int bpp = 0;
  int flags = 0;

  if( SDL_Init( SDL_INIT_VIDEO ) < 0 )
  {
    fprintf( stderr, "[display] video initialization failed: %s\n", SDL_GetError() );
    free(d);
    return 0;
  }

  info = SDL_GetVideoInfo( );

  if( !info )
  {
    fprintf( stderr, "[display] video query failed: %s\n", SDL_GetError() );
    free(d);
    return 0;
  }
  bpp = 32;

  SDL_GL_SetAttribute( SDL_GL_RED_SIZE, 8 );
  SDL_GL_SetAttribute( SDL_GL_GREEN_SIZE, 8 );
  SDL_GL_SetAttribute( SDL_GL_BLUE_SIZE, 8 );
  SDL_GL_SetAttribute( SDL_GL_DEPTH_SIZE, 16 );
  SDL_GL_SetAttribute( SDL_GL_DOUBLEBUFFER, 1 );

  flags = SDL_OPENGL;// | SDL_FULLSCREEN;

  if( SDL_SetVideoMode( width, height, bpp, flags ) == 0 )
  {
    fprintf( stderr, "[display] video mode set failed: %s\n", SDL_GetError( ) );
    free(d);
    return 0;
  }
  SDL_WM_SetCaption(title, NULL);

  atexit(&SDL_Quit);

#if 0
  GLenum err = glewInit();
  if (err != GLEW_OK)
  {
    fprintf(stderr, "[display] error: %s\n", glewGetErrorString(err));
    return 0;
  }
  printf("[display] Using glew %s\n", glewGetString(GLEW_VERSION));

  glClearColor (0.0, 0.0, 0.0, 0.0);
  glShadeModel(GL_FLAT);
  glEnable(GL_TEXTURE_2D);

  GLuint buf_size = width * height * 3;
  GLuint id;

  glGenBuffers(1, &id);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, 0);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

  glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, id);
  glBufferData(GL_PIXEL_UNPACK_BUFFER_ARB, buf_size, 0, GL_STREAM_DRAW);
#endif

#if 1//def TEXTURE
  // printf("blitting using tex image\n");
  GLuint texID = 0;
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_BLEND);
  glColor3f(1.0f, 1.0f, 1.0f);
  glGenTextures(1, &texID);
  glBindTexture(GL_TEXTURE_2D, texID);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
  //glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, WIDTH, HEIGHT, 0, GL_RGB, GL_INT, NULL);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, d->width, d->height, 0, GL_RGBA, GL_FLOAT, NULL);
  
  glEnable(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D, texID);

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();  glLoadIdentity();
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
#endif

  printf("[display] keyboard control:");
#ifdef GUI_QWERTZ
  printf(" (qwertz)\n");
  printf(" [ye]  camera speed -/+\n");
#else
  printf(" (dvorak)\n");
  printf(" [;.]  camera speed -/+\n");
#endif
  printf(" [12]  exposure time\n");
  printf(" [34]  f-stop\n");
  printf(" [56]  focal len\n");
  printf(" [78]  iso value\n");
  printf(" [c]   save camera position\n");
  printf(" [l]   load camera position\n");
  printf(" [p]   print screenshot\n");
  printf(" [h]   toggle display\n");
  printf(" quake movement, use left mouse button to turn\n");

  return d;
}

void display_close(display_t *d)
{	
  free(d);
}

void handleEvent(const SDL_Event *event, display_t *d)
{
  switch (event->type)
  {
    case SDL_KEYDOWN:
      {
      const SDLKey keysym = event->key.keysym.sym;
      if(keysym == SDLK_ESCAPE)
      {
        d->isShuttingDown = 1;
        if(d->onClose) d->onClose();
      }
      else if(keysym < keyMapSize && d->onKeyDown) d->onKeyDown(normalKeys[keysym]);
      }
      break;
    case SDL_KEYUP:
      {
      const SDLKey keysym = event->key.keysym.sym;
      if(keysym < keyMapSize && d->onKeyUp) d->onKeyUp(normalKeys[keysym]);
      }
      break;
    case SDL_QUIT:
      d->isShuttingDown = 1;
      if(d->onClose) d->onClose();
      break;
    case SDL_MOUSEBUTTONDOWN:
    case SDL_MOUSEBUTTONUP:
      {
        mouse_t mouse;
        mouse.x = event->button.x;
        mouse.y = event->button.y;
        mouse.buttons.left   = event->button.button & SDL_BUTTON_LEFT;
        mouse.buttons.middle = event->button.button & SDL_BUTTON_MIDDLE;
        mouse.buttons.right  = event->button.button & SDL_BUTTON_RIGHT;
        if (event->type == SDL_MOUSEBUTTONDOWN)
        {
          if (d->onMouseButtonDown) d->onMouseButtonDown(mouse);
        }
        else if (d->onMouseButtonUp) d->onMouseButtonUp(mouse);
        break;
      }
    case SDL_MOUSEMOTION:
      {
        mouse_t mouse;
        mouse.x = event->motion.x;
        mouse.y = event->motion.y;
        mouse.buttons.left   = event->button.button & SDL_BUTTON_LEFT;
        mouse.buttons.middle = event->button.button & SDL_BUTTON_MIDDLE;
        mouse.buttons.right  = event->button.button & SDL_BUTTON_RIGHT;
        if (d->onMouseMove) d->onMouseMove(mouse);
        break;
      }
  }
}

void display_pump_events(display_t *d)
{
  SDL_Event event;
  while(SDL_PollEvent(&event)) handleEvent(&event, d);
}

void convert_3(unsigned char* restrict bbuf, const float* restrict fbuf, int size)
{
  const int count = (((size)>>2)&0xFFFFFFFC);
  const __m128 c256 = _mm_set1_ps(256.0f);
  const __m128 c1 = _mm_set1_ps(1.0f);
  const __m128* f = (__m128*)fbuf;
  //__m128i* c = (__m128i*)bbuf;
// TODO: port to thread pool:
// #pragma omp parallel for default(none) schedule(static) shared(bbuf, fbuf, f)
  for (int j=0; j < (count>>2); j++) { // converts 16 floats at once
    const int i = j<<2;
    __m128i *c = ((__m128i*)bbuf) + j;
    //_mm_prefetch((void*)(f+i+32), _MM_HINT_T0); // two cache lines ahead
    // gamma = 2.0
    // __m128i f0 =       _mm_cvttps_epi32(_mm_mul_ps(_mm_mul_ps(f[i+0], f[i+0]), c256));
    // const __m128i f1 = _mm_cvttps_epi32(_mm_mul_ps(_mm_mul_ps(f[i+1], f[i+1]), c256));
    // __m128i f2 =       _mm_cvttps_epi32(_mm_mul_ps(_mm_mul_ps(f[i+2], f[i+2]), c256));
    // const __m128i f3 = _mm_cvttps_epi32(_mm_mul_ps(_mm_mul_ps(f[i+3], f[i+3]), c256));
#ifdef MACOSX
    __m128i f0 =       _mm_cvttps_epi32(_mm_max_ps(_mm_setzero_ps(), _mm_mul_ps(_mm_mul_ps(f[i+0], _mm_rcp_ps(_mm_add_ps(f[i+0], c1))), c256)));
    const __m128i f1 = _mm_cvttps_epi32(_mm_max_ps(_mm_setzero_ps(), _mm_mul_ps(_mm_mul_ps(f[i+1], _mm_rcp_ps(_mm_add_ps(f[i+1], c1))), c256)));
    __m128i f2 =       _mm_cvttps_epi32(_mm_max_ps(_mm_setzero_ps(), _mm_mul_ps(_mm_mul_ps(f[i+2], _mm_rcp_ps(_mm_add_ps(f[i+2], c1))), c256)));
    const __m128i f3 = _mm_cvttps_epi32(_mm_max_ps(_mm_setzero_ps(), _mm_mul_ps(_mm_mul_ps(f[i+3], _mm_rcp_ps(_mm_add_ps(f[i+3], c1))), c256)));
    f0 = _mm_packs_epi32(_mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(f0), _mm_castsi128_ps(f0), _MM_SHUFFLE(0,1,2,3))), _mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(f1), _mm_castsi128_ps(f1), _MM_SHUFFLE(0,1,2,3))));
    f2 = _mm_packs_epi32(_mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(f2), _mm_castsi128_ps(f2), _MM_SHUFFLE(0,1,2,3))), _mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(f3), _mm_castsi128_ps(f3), _MM_SHUFFLE(0,1,2,3))));
#else
    // tone mapping:
    __m128i f0 =       _mm_cvttps_epi32(_mm_mul_ps(_mm_mul_ps(f[i+0], _mm_rcp_ps(_mm_add_ps(f[i+0], c1))), c256));
    const __m128i f1 = _mm_cvttps_epi32(_mm_mul_ps(_mm_mul_ps(f[i+1], _mm_rcp_ps(_mm_add_ps(f[i+1], c1))), c256));
    __m128i f2 =       _mm_cvttps_epi32(_mm_mul_ps(_mm_mul_ps(f[i+2], _mm_rcp_ps(_mm_add_ps(f[i+2], c1))), c256));
    const __m128i f3 = _mm_cvttps_epi32(_mm_mul_ps(_mm_mul_ps(f[i+3], _mm_rcp_ps(_mm_add_ps(f[i+3], c1))), c256));
    f0 = _mm_packs_epi32(f0, f1);
    f2 = _mm_packs_epi32(f2, f3);
#endif
    __m128i result = _mm_packus_epi16(f0, f2); // saturates at 255
    //_mm_stream_si128(c, result);
    *c = result;
  }
  for(int i=count*4; i < size; i++) bbuf[i] = fminf(1.f,(fbuf[i]))*255.0f; // convert the rest
  _mm_mfence();
}

int display_update(display_t *d, float * pixels)
{
  if (d->isShuttingDown)
  {
    display_close(d);
    return 0;
  }
#if 1
  // render message:
  int px = d->msg_x;
  for (int pos = 0; pos < d->msg_len; pos++)
  {
    int charPos = (d->msg[pos] - 32)*9*2;
    for (int x = 0; x < 9*2; x++)
    {
      if (px >= d->width) goto render_message_out;
      unsigned char cLine = font9x16[charPos+x];			
      for (int i = 0; i < 8; i++)
      {
        int y = d->msg_y - (15 - (i + (x&1)*8));
        if ((y >= d->height) || (y < 0)) goto render_message_out;
        if (cLine & (1<<(7-i)))
          for(int k=1;k<4;k++) pixels[(px + y*d->width)*4+k] = .7;
        else
          for(int k=1;k<4;k++) pixels[(px + y*d->width)*4+k] = CLAMP(pixels[(px + y*d->width)*4+k], 0.0f, 1.0f)*0.3f;
      }
      if (x&1) px++;
    }
  }
render_message_out:
#endif
#if 1//def TEXTURE
  glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, d->width, d->height, GL_RGBA, GL_FLOAT, pixels+1); // alles ist scheisse.
  glBegin(GL_QUADS);
  glTexCoord2f(0.0, 1.0); glVertex3f(-1.0, -1.0, 0.0);
  glTexCoord2f(0.0, 0.0); glVertex3f(-1.0, 1.0, 0.0);
  glTexCoord2f(1.0, 0.0); glVertex3f(1.0, 1.0, 0.0);
  glTexCoord2f(1.0, 1.0); glVertex3f(1.0, -1.0, 0.0);
  glEnd();
#endif

#if 0//def FBO
  mem = (GLubyte *)glMapBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, GL_WRITE_ONLY);
  //if(!mem) exit(1);//continue;
  for (unsigned int y = 0; y < d->height; ++y)
    for (unsigned int x = 0; x < d->width ; ++x) {
      mem[y*d->width*3 + x*3] = y % 255;
      mem[y*d->width*3 + x*3+1] = sweep%255;
      mem[y*d->width*3 + x*3+2] = x % 255;
    }
  glUnmapBuffer(GL_PIXEL_UNPACK_BUFFER_ARB);
  glTexSubImage2D(GL_TEXTURE_2D, 0, 0,0, d->width, d->height, GL_RGB, GL_UNSIGNED_BYTE, 0);

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glBegin(GL_QUADS);
  glTexCoord2f(0.0, 0.0); glVertex3f(-1.0, -1.0, 0.0);
  glTexCoord2f(0.0, 1.0); glVertex3f(-1.0, 1.0, 0.0);
  glTexCoord2f(1.0, 1.0); glVertex3f(1.0, 1.0, 0.0);
  glTexCoord2f(1.0, 0.0); glVertex3f(1.0, -1.0, 0.0);
  glEnd();
  glFlush();
#endif

  SDL_GL_SwapBuffers();
  return 1;
}

int display_update_rgba(display_t *d, const unsigned int * rgba)
{
  if (d->isShuttingDown)
  {
    display_close(d);
    return 0;
  }

  return 1;
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
