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

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <dlfcn.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <assert.h>

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

static keycode_t normalKeys[keyMapSize];
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
	for (int i = 0; i < keyMapSize; ++i) {
		normalKeys[i] = KeyUndefined;
		// functionKeys[i] = KeyUndefined;
	}

	normalKeys[GLFW_KEY_SPACE] = KeySpace;
	normalKeys[GLFW_KEY_COMMA] = KeyComma;
	normalKeys[GLFW_KEY_PERIOD] = KeyPeriod;
	normalKeys[GLFW_KEY_SLASH] = KeySlash;
	normalKeys[GLFW_KEY_0] = KeyZero;
	normalKeys[GLFW_KEY_1] = KeyOne;
	normalKeys[GLFW_KEY_2] = KeyTwo;
	normalKeys[GLFW_KEY_3] = KeyThree;
	normalKeys[GLFW_KEY_4] = KeyFour;
	normalKeys[GLFW_KEY_5] = KeyFive;
	normalKeys[GLFW_KEY_6] = KeySix;
	normalKeys[GLFW_KEY_7] = KeySeven;
	normalKeys[GLFW_KEY_8] = KeyEight;
	normalKeys[GLFW_KEY_9] = KeyNine;
	normalKeys[GLFW_KEY_SEMICOLON] = KeySemiColon;
	normalKeys[GLFW_KEY_EQUAL] = KeyEquals;
	normalKeys[GLFW_KEY_A] = KeyA;
	normalKeys[GLFW_KEY_B] = KeyB;
	normalKeys[GLFW_KEY_C] = KeyC;
	normalKeys[GLFW_KEY_D] = KeyD;
	normalKeys[GLFW_KEY_E] = KeyE;
	normalKeys[GLFW_KEY_F] = KeyF;
	normalKeys[GLFW_KEY_G] = KeyG;
	normalKeys[GLFW_KEY_H] = KeyH;
	normalKeys[GLFW_KEY_I] = KeyI;
	normalKeys[GLFW_KEY_J] = KeyJ;
	normalKeys[GLFW_KEY_K] = KeyK;
	normalKeys[GLFW_KEY_L] = KeyL;
	normalKeys[GLFW_KEY_M] = KeyM;
	normalKeys[GLFW_KEY_N] = KeyN;
	normalKeys[GLFW_KEY_O] = KeyO;
	normalKeys[GLFW_KEY_P] = KeyP;
	normalKeys[GLFW_KEY_Q] = KeyQ;
	normalKeys[GLFW_KEY_R] = KeyR;
	normalKeys[GLFW_KEY_S] = KeyS;
	normalKeys[GLFW_KEY_T] = KeyT;
	normalKeys[GLFW_KEY_U] = KeyU;
	normalKeys[GLFW_KEY_V] = KeyV;
	normalKeys[GLFW_KEY_W] = KeyW;
	normalKeys[GLFW_KEY_X] = KeyX;
	normalKeys[GLFW_KEY_Y] = KeyY;
	normalKeys[GLFW_KEY_Z] = KeyZ;
	normalKeys[GLFW_KEY_LEFT_BRACKET] = KeyOpenBracket;
	normalKeys[GLFW_KEY_BACKSLASH] = KeyBackSlash;
	normalKeys[GLFW_KEY_RIGHT_BRACKET] = KeyCloseBracket;
	return 1;
}

static void
glfw_error_callback(int error, const char *description)
{
	fprintf(stderr, "glfw error [%d]: %s\n", error, description);
}

static void
glfw_key_callback(GLFWwindow *win, int key, int scancode, int action, int mods)
{
	display_t *d = glfwGetWindowUserPointer(win);
	switch(key) {
	case GLFW_KEY_ESCAPE:
		glfwSetWindowShouldClose(win, GL_TRUE);
		d->isShuttingDown = 1;
		if(d->onClose)
			d->onClose();
		break;
	default:
		if(action == GLFW_PRESS)
    {
			if(key < keyMapSize && d->onKeyDown)
				d->onKeyDown(normalKeys[key]);
		}
		else if(action == GLFW_RELEASE)
    {
			if(key < keyMapSize && d->onKeyUp)
				d->onKeyUp(normalKeys[key]);
		}
		break;
	}
}

static void
glfw_button_callback(GLFWwindow *win, int button, int action, int mods)
{
	display_t *d = glfwGetWindowUserPointer(win);
	double x, y;
	glfwGetCursorPos(win, &x, &y);
	mouse_t mouse;
	mouse.x = x;
	mouse.y = y;
	mouse.buttons.left   = glfwGetMouseButton(win, 0);
	mouse.buttons.middle = glfwGetMouseButton(win, 1);
	mouse.buttons.right  = glfwGetMouseButton(win, 2);
	if(action == GLFW_PRESS) {
		if (d->onMouseButtonDown) d->onMouseButtonDown(mouse);
	}
	else {
		if (d->onMouseButtonUp) d->onMouseButtonUp(mouse);
	}
}

static void
glfw_motion_callback(GLFWwindow *win, double x, double y)
{
	display_t *d = glfwGetWindowUserPointer(win);

	mouse_t mouse;
	mouse.x = x;
	mouse.y = y;
	mouse.buttons.left   = glfwGetMouseButton(win, 0);
	mouse.buttons.middle = glfwGetMouseButton(win, 1);
	mouse.buttons.right  = glfwGetMouseButton(win, 1);
	if (d->onMouseMove)
		d->onMouseMove(mouse);
}

static unsigned int
compile_src(const char *src, int type)
{
	GLuint shader_obj = glCreateShader(type);

	glShaderSource(shader_obj, 1, &src, NULL);
	glCompileShader(shader_obj);

	int success;
	glGetShaderiv(shader_obj, GL_COMPILE_STATUS, &success);
	if(!success) {
		char infolog[2048];
		glGetShaderInfoLog(shader_obj, sizeof infolog, NULL, infolog);
		fprintf(stderr, "error compiling shader: %s\n", infolog);
	}

	return shader_obj;
}

static unsigned int
compile_shader(const char *src_vertex, const char *src_fragment)
{
	GLuint prog = glCreateProgram();

	unsigned int vert = 0, frag = 0;

	glAttachShader(prog, vert = compile_src(src_vertex, GL_VERTEX_SHADER));
	glAttachShader(prog, frag = compile_src(src_fragment, GL_FRAGMENT_SHADER));
	glLinkProgram(prog);

	int success;
	glGetProgramiv(prog, GL_LINK_STATUS, &success);
	if(!success) {
		char infolog[2048];
		glGetProgramInfoLog(prog, sizeof infolog, NULL, infolog);
		fprintf(stderr, "error linking shaderprogram: %s\n", infolog);
	}

	glDeleteShader(vert);
	glDeleteShader(frag);

	return prog;
}

static void
gl_error_callback(GLenum source, GLenum type, GLuint id, GLenum severity,
		GLsizei length, const char *msg, GLvoid *param)
{
	const char *ssource;
	switch(source) {
	case 0x8246: ssource = "API            "; break;
	case 0x8247: ssource = "WINDOW_SYSTEM  "; break;
	case 0x8248: ssource = "SHADER_COMPILER"; break;
	case 0x8249: ssource = "THIRD_PARTY    "; break;
	case 0x824A: ssource = "APPLICATION    "; break;
	case 0x824B: ssource = "OTHER          "; break;
	default:     ssource = "INVALID        "; break;
	}

	const char *stype;

	switch(type) {
	case 0x824C: stype = "ERROR              "; break;
	case 0x824D: stype = "DEPRECATED_BEHAVIOR"; break;
	case 0x824E: stype = "UNDEFINED_BEHAVIOR "; break;
	case 0x824F: stype = "PORTABILITY        "; break;
	case 0x8250: stype = "PERFORMANCE        "; break;
	case 0x8251: stype = "OTHER              "; break;
	default:     stype = "INVALID            "; break;
	}

	const char *sseverity;

	switch(severity) {
	case 0x9146: sseverity = "HIGH   "; break;
	case 0x9147: sseverity = "MEDIUM "; break;
	case 0x9148: sseverity = "LOW    "; break;
	default:     sseverity = "INVALID"; break;
	}

	fprintf(stderr, "%s %s %s: %s\n", sseverity, stype, ssource, msg);
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

	glfwSetErrorCallback(glfw_error_callback);
	if(!glfwInit()) {
		fprintf(stderr, "error initializing glfw\n");
		free(d);
		return NULL;
	}

	glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
	// glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	// glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3); // seems to be the newest the quadro 5800 supports :(
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_RED_BITS, 10);
  glfwWindowHint(GLFW_GREEN_BITS, 10);
  glfwWindowHint(GLFW_BLUE_BITS, 10);
  // glfwWindowHint(GLFW_ALPHA_BITS, 2);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	if(!(d->window = glfwCreateWindow(width, height, title, NULL, NULL))) {
		fprintf(stderr, "error creating window\n");
		glfwTerminate();
		free(d);
		return NULL;
	}

	glfwMakeContextCurrent(d->window);
	glfwSetKeyCallback(d->window, glfw_key_callback);
	glfwSetCursorPosCallback(d->window, glfw_motion_callback);
	glfwSetMouseButtonCallback(d->window, glfw_button_callback);

	glewExperimental = GL_TRUE;
	GLenum err = glewInit();
	if(err != GLEW_OK) {
		fprintf(stderr,  "Error: %s\n", glewGetErrorString(err));
		return NULL;
	}
	while(glGetError() != GL_NO_ERROR);
	glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
	glEnable(GL_DEBUG_OUTPUT);

	glDebugMessageCallback((GLDEBUGPROC) gl_error_callback, NULL);
	//glDebugMessageControl(GL_DONT_CARE, GL_DONT_CARE, GL_DONT_CARE, 0, NULL, GL_TRUE);
	GLuint mask_ids[] = {
		131076, // Generic Vertex Attrib foo
		131185, // memory type for buffer operations
	};
	glDebugMessageControl(GL_DEBUG_SOURCE_API, GL_DEBUG_TYPE_OTHER, GL_DONT_CARE,
			sizeof(mask_ids) / sizeof(mask_ids[0]), mask_ids, GL_FALSE);

	d->program_draw_texture = compile_shader(
"#version 330\n"
"out vec2 tex_coord;\n"
"void\n"
"main()\n"
"{\n"
"	vec2 p;"
"	if(gl_VertexID == 0)\n"
"		p = vec2(-1,  1);\n"
"	else if(gl_VertexID == 2)\n"
"		p = vec2(-1, -3);\n"
"	else\n"
"		p = vec2( 3,  1);\n"
"	tex_coord = p.xy * 0.5 + 0.5;\n"
"	gl_Position = vec4(p, 0.0, 1.0);\n"
"}\n",

"#version 330\n"
"in vec2 tex_coord;\n"
"uniform sampler2D framebuffer;\n"
"out vec4 frag_color;\n"
"void\n"
"main()\n"
"{\n"
"	frag_color = texture(framebuffer, vec2(tex_coord.x, 1.0 - tex_coord.y));\n"
"}\n"
);
	glGenVertexArrays(1, &d->vao_empty);
	glGenTextures(1, &d->tex_frame_buffer);
	glBindTexture(GL_TEXTURE_2D, d->tex_frame_buffer);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, width, height, 0,
		GL_RGBA, GL_FLOAT, NULL);
	glGenerateMipmap(GL_TEXTURE_2D);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

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
	glfwDestroyWindow(d->window);
	glfwTerminate();

	free(d);
}

void display_pump_events(display_t *d)
{
	glfwSetWindowUserPointer(d->window, d);
	glfwPollEvents();

	if(glfwWindowShouldClose(d->window)) {
		d->isShuttingDown = 1;
		if(d->onClose)
			d->onClose();
	}
}

int display_update(display_t *d, float * pixels)
{
	if (d->isShuttingDown)
	{
		display_close(d);
		return 0;
	}

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

	glClearColor(0, 1, 0, 1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glUseProgram(d->program_draw_texture);

	glBindTexture(GL_TEXTURE_2D, d->tex_frame_buffer);
	glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, d->width, d->height, GL_RGBA, GL_FLOAT, pixels + 1);

	glBindVertexArray(d->vao_empty);
	glDrawArrays(GL_TRIANGLE_STRIP, 0, 3);

	glfwSwapBuffers(d->window);
	return 1;
}

void display_print(display_t *d, const int px, const int py, const char *msg, ...)
{
	va_list ap;
	va_start(ap, msg);
	vsnprintf(d->msg, 255, msg, ap);
	va_end(ap);
	d->msg_len = strlen(d->msg);
	d->msg_x = px;
	d->msg_y = py + 15;
}
