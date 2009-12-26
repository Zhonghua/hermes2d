// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#ifndef NOGLUT

#include <GL/freeglut.h>
#ifndef WIN32
# include <sys/time.h>
#endif

#include "../common.h"
#include "view_support.h"
#include "view.h"
#include "../solution.h"

///////////////// static variables /////////////////
int View::screenshot_no = 1;
View::ViewMonitor View::view_sync;

///////////////// monitor /////////////////
View::ViewMonitor::ViewMonitor()
{
  pthread_mutex_init(&mutex, NULL);
  pthread_cond_init(&cond_keypress, NULL);
  pthread_cond_init(&cond_close, NULL);
  pthread_cond_init(&cond_drawing_finished, NULL);
}

View::ViewMonitor::~ViewMonitor()
{
  pthread_mutex_destroy(&mutex);
  pthread_cond_destroy(&cond_keypress);
  pthread_cond_destroy(&cond_close);
  pthread_cond_destroy(&cond_drawing_finished);
}

///////////////// methods /////////////////
View::View(const char* title, int x, int y, int width, int height)
  : gl_pallete_tex_id(0)
{
  this->title = title;
  window_x = x;
  window_y = y;
  window_width = width;
  window_height = height;
  window_id = -1;
  jitter_x = jitter_y = 0.0;
  dragging = scaling = false;
  hq_frame = false;
  frame_ready = false;
  range_auto = true;
  range_min = 0;
  range_max = 1;
  pal_type = 0;
  pal_steps = 50;
  pal_filter = GL_NEAREST;
  margin = 15;
  b_scale = true;
  b_help = false;
  scale_focused = scale_dragging = false;
  pos_horz = pos_vert = 0;
  scale_x = scale_y = labels_width = 0;
  scale_width = 16;
  scale_height = 320;
  scale_numticks = 9;
  strcpy(scale_fmt, "%.3g");
  scale_fixed_width = -1;
  want_screenshot = false;
}


View::~View()
{
  if (window_id >= 0)
  {
    verbose("View is being destroyed; closing window #%d.", window_id);
    close();
    if (!num_windows) finish_glut_main_loop(false);
  }
}


int view_create_body(void* param)
{
  View* instance = (View*) param;

  // create the window
  glutInitWindowPosition(instance->window_x, instance->window_y);
  glutInitWindowSize(instance->window_width, instance->window_height);
  instance->window_id = glutCreateWindow(instance->title.c_str());
  num_windows++;

  // establish a mapping between the window id and this instance
  if (instance->window_id >= (int)wnd_instance.size())
    wnd_instance.resize(instance->window_id + 10, NULL);
  wnd_instance[instance->window_id] = instance;

  // register callbacks
  glutDisplayFunc(on_display_stub);
  glutReshapeFunc(on_reshape_stub);
  glutMotionFunc(on_mouse_move_stub);
  glutPassiveMotionFunc(on_mouse_move_stub);
  glutMouseFunc(on_mouse_click_stub);
  glutKeyboardFunc(on_key_down_stub);
  glutSpecialFunc(on_special_key_stub);
  glutEntryFunc(on_entry_stub);
  glutCloseFunc(on_close_stub);

  instance->on_create();
  return instance->window_id;
}

void View::init_output()
{
  glut_init();
}

int View::create()
{
  //reset timing
  memset(rendering_frames, 0, FPS_FRAME_SIZE * sizeof(double));
  rendering_frames_top = 0;

  if (window_id >= 0)
    safe_post_redisplay();
  else
    cross_thread_call(view_create_body, this);
  return window_id;
}


void View::close()
{
  if (window_id >= 0)
  {
    glutDestroyWindow(window_id);
    if (!pthread_equal(thread, pthread_self())) wait_for_close();
    window_id = -1;
  }
}


void View::wait(const char* text)
{
  if (text != NULL) printf("%s\n", text);
  finish_glut_main_loop(false);
}


void View::on_create()
{
  create_gl_palette();
  set_palette_filter(pal_filter == GL_LINEAR);
}


void View::on_close()
{
  verbose("Window #%d closed.", glutGetWindow());
  view_sync.enter();
  view_sync.signal_close();
  view_sync.signal_keypress();
  view_sync.signal_drawing_finished();
  view_sync.leave();
}


void View::clear_background()
{
  glClearColor(1.0, 1.0, 1.0, 0.0);
  glClear(GL_COLOR_BUFFER_BIT);
}


void View::pre_display()
{
  //info("display: lock");
  view_sync.enter();

  //begin time measuring
  double time_start = get_tick_count();

//  if (!hq_frame)
//  {
    clear_background();
    on_display();
//  }
//  else 
//  {
//    display_antialiased();
//    hq_frame = false;
//  }
  
  if (b_help) draw_help();
  else if (b_scale) scale_dispatch();

  //draw current rendring time
  draw_fps();

  //wait to finish
  glFinish();

  //calculate statistics
  rendering_frames[rendering_frames_top] = get_tick_count() - time_start;
  rendering_frames_top = (rendering_frames_top + 1) % FPS_FRAME_SIZE;

  if (want_screenshot)
  {
    glReadBuffer(GL_BACK_LEFT);
    save_screenshot_internal(screenshot_filename.c_str());
    want_screenshot = false;
  }

  glutSwapBuffers();

  //frame synchronization
  frame_ready = true;
  view_sync.signal_drawing_finished();
  view_sync.leave();
}


static float jitter16[16][2] =
{
  { 0.4375, 0.4375 }, { 0.1875, 0.5625 },
  { 0.9375, 1.1875 }, { 0.4375, 0.9375-1 },
  { 0.6875, 0.5625 }, { 0.1875, 0.0625 },
  { 0.6875, 0.3125 }, { 0.1875, 0.3125 },
  { 0.4375, 0.1875 }, { 0.9375-1, 0.4375 },
  { 0.6875, 0.8125 }, { 0.4375, 0.6875 },
  { 0.6875, 0.0625 }, { 0.9375, 0.9375 },
  { 1.1875, 0.8125 }, { 0.9375, 0.6875 }
};


void View::display_antialiased()
{
  glClear(GL_ACCUM_BUFFER_BIT);
  for (int i = 0; i < 16; i++)
  {
    jitter_x = jitter16[i][0];
    jitter_y = jitter16[i][1];
    set_ortho_projection();
    clear_background();
    on_display();
    glAccum(GL_ACCUM, 1.0 / 16);
  }
  glAccum(GL_RETURN, 1.0);
  jitter_x = jitter_y = 0.0;
}


void View::set_ortho_projection(bool no_jitter)
{
  double jx = no_jitter ? 0.0 : jitter_x;
  double jy = no_jitter ? 0.0 : jitter_y;

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(jx, window_width + jx, window_height-1 + jy, -1 + jy, -10, 10);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

void View::set_3d_projection(int fov, double znear, double zfar)
{
  double right = znear * tan((double) fov / 2.0 / 180.0 * M_PI);
  double top = (double) window_height / window_width * right;
  double left = -right;
  double bottom = -top;
	double offsx = (right - left) / window_width * jitter_x;
	double offsy = (top - bottom) / window_height * jitter_y;

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glFrustum(left - offsx, right - offsx, bottom - offsy, top - offsy, znear, zfar);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

void View::draw_fps()
{
  //calculate FPS
  double frame_time_sum = 0;
  for(int i = 0; i < FPS_FRAME_SIZE; i++)
    frame_time_sum += rendering_frames[i];

  //prepare text
  unsigned char buffer[128];
  sprintf((char*)buffer, "avg. frame: %.1f ms", (float)(frame_time_sum / FPS_FRAME_SIZE));

  //prepare environment
  void* font = GLUT_BITMAP_HELVETICA_10;
  int width_px = glutBitmapLength(font, buffer);
  int height_px = glutBitmapHeight(font);
  int edge_thickness = 2;
  set_ortho_projection(false);

  //render background
  glEnable(GL_BLEND);
  glBlendFunc(GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
  glBegin(GL_QUADS);
  glColor4f(1.0f, 1.0f, 1.0f, 0.5f);
  glVertex2i(window_width - (width_px + 2*edge_thickness), 0);
  glVertex2i(window_width, 0);
  glVertex2i(window_width, height_px + 2*edge_thickness);
  glVertex2i(window_width - (width_px + 2*edge_thickness), height_px + 2*edge_thickness);
  glEnd();

  //render text
  glDisable(GL_BLEND);
  glColor3f(1.0f, 0.0f, 0.0f);
  glRasterPos2i(window_width - (width_px + edge_thickness), edge_thickness + height_px);
  glutBitmapString(font, buffer);
}

void View::on_reshape(int width, int height)
{
  window_width = width;
  window_height = height;
  update_layout();
  glViewport(0, 0, width, height);

  /*printf("winx=%d, winy=%d, ww=%d, wh=%d, width=%d, height=%d\n", glutGet(GLUT_WINDOW_X), glutGet(GLUT_WINDOW_Y),
         glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT), width, height);*/
}


void View::update_scale()
{
  scale = pow(1.005, log_scale);
}

void View::update_log_scale()
{
  log_scale = log(scale) / log(1.005);
}


void View::on_mouse_move(int x, int y)
{
  if (dragging)
  {
    trans_x += (x - mouse_x);
    trans_y += (mouse_y - y);
    post_redisplay();
  }
  else if (scaling)
  {
    log_scale += (mouse_y - y);
    update_scale();
    trans_x = scx - objx * scale - center_x;
    trans_y = center_y - objy * scale - scy;
    //on_zoom(get_zoom());
    post_redisplay();
  }
  else if (scale_dragging)
  {
    int oldv = pos_vert, oldh = pos_horz;
    pos_horz = (x > window_width/2);
    pos_vert = (y < window_height/2);
    if (pos_horz != oldh || pos_vert != oldv) { update_layout(); post_redisplay(); }
  }
  else
  {
    bool oldf = scale_focused;
    scale_focused = (x >= scale_x && x <= scale_x + scale_width &&
                     y >= scale_y && y <= scale_y + scale_height);
    if (oldf != scale_focused) post_redisplay();
  }
  mouse_x = x;
  mouse_y = y;
}


void View::on_left_mouse_down(int x, int y)
{
  if (scale_focused)
    scale_dragging = true;
  else
    dragging = true;
  scaling = false;
  mouse_x = x;
  mouse_y = y;
}


void View::on_left_mouse_up(int x, int y)
{
  scaling = dragging = scale_dragging = false;
  on_mouse_move(x, y);
}


void View::on_right_mouse_down(int x, int y)
{
  scaling = true;
  dragging = false;
  scx = x;
  scy = y;
  objx = (x - center_x - trans_x) / scale;
  objy = (center_y - y - trans_y) / scale;
  mouse_x = x;
  mouse_y = y;
}


void View::on_right_mouse_up(int x, int y)
{
  scaling = dragging = false;
}


void View::on_key_down(unsigned char key, int x, int y)
{
  const char *file_name = NULL;

  switch (key)
  {
    case 'h':
    {
      hq_frame = true;
      post_redisplay();
      break;
    }

    case 27:
    case 'q':
    {
      close();
      break;
    }

    case 's':
    {
      const char *file_name = get_screenshot_file_name();
      glReadBuffer(GL_FRONT_LEFT);
      save_screenshot_internal(file_name);
      break;
    }

    case 'p':
    {
      pal_type++;
      if (pal_type > 3) pal_type = 0; 
      create_gl_palette();
      post_redisplay();
      break;
    }

    default:
      view_sync.enter();
      view_sync.signal_keypress();
      view_sync.leave();
      break;
  }
}


void View::on_special_key(int key, int x, int y)
{
  switch (key)
  {
    case GLUT_KEY_F1:
      b_help = !b_help;
      post_redisplay();
      break;
  }
}


void View::wait_for_keypress()
{
  view_sync.enter();
  if (window_id >= 0)
    view_sync.wait_keypress();
  view_sync.leave();
}

void View::wait_for_close()
{
  view_sync.enter();
  if (window_id >= 0)
    view_sync.wait_close();
  view_sync.leave();
}

void View::wait_for_draw()
{
  if (!pthread_equal(thread, pthread_self()))
  {
    view_sync.enter();
    if (window_id >= 0 && !frame_ready)
      view_sync.wait_drawing_fisnihed();
    view_sync.leave();
  }
  //old code
  //if (window_id < 0) return;
  //if (pthread_equal(thread, pthread_self())) return; //?WTF: if my thread is not my thread then lock and wait? This means that the function can be called by the thread that is drawing?
  ////info("wait: lock");
  //pthread_mutex_lock(&wait_draw_mutex);
  //if (!frame_ready) { //info("wait: wait");
  //  pthread_cond_wait(&wait_draw_cv, &wait_draw_mutex); }
  ////info("wait: unlock");
  //pthread_mutex_unlock(&wait_draw_mutex);
}


void View::post_redisplay()
{
  debug_assert(pthread_equal(thread, pthread_self()) != 0, "E 'View::post_redisplay' accessed from other than drawing thread.\n");
  if (window_id < 0) return; //set by this thread, other threads just reads
  //info("post: lock");
  //pthread_mutex_lock(&wait_draw_mutex);
  glutPostWindowRedisplay(window_id);
  frame_ready = false;
  //info("post");
  //info("post: unlock");
  //pthread_mutex_unlock(&wait_draw_mutex);
}


static int safe_post_redisplay_body(void* param)
{
  glutPostWindowRedisplay((int) (long) param);
  return 0;
}

void View::safe_post_redisplay()
{
  if (!pthread_equal(thread, pthread_self()))
  {
    if (window_id >= 0)
      cross_thread_call(safe_post_redisplay_body, (void*) window_id);
  }
  else {
    post_redisplay();
  }
}


double View::get_tick_count()
{
#ifdef WIN32
  LARGE_INTEGER freq, ticks;
  QueryPerformanceFrequency(&freq);
  QueryPerformanceCounter(&ticks);
  return (1000.0 * (double)ticks.QuadPart) / (double)freq.QuadPart;
  //return (double) clock() * (1000.0 / (double) CLOCKS_PER_SEC);
#else
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (double) tv.tv_sec * 1000 + (double) tv.tv_usec / 1000;
#endif
}


int view_set_title_body(void* param)
{
  View *self = (View*) param;
  self->set_title_internal(self->title.c_str());

  return 0;
}

void View::set_title(const char* title)
{
  this->title = title;
  if (window_id >= 0)
    cross_thread_call(view_set_title_body, this);
}


void View::set_title_internal(const char* text)
{
  glutSetWindow(window_id);
  glutSetWindowTitle(text);
}


//// palette ///////////////////////////////////////////////////////////////////////////////////////

#include "view_data.cpp"

const float* View::get_palette_color(double x)
{
  static float color[3];

  if (pal_type == 0)
  {
    if (x < 0.0) x = 0.0;
    else if (x > 1.0) x = 1.0;
    x *= num_pal_entries;
    int n = (int) x;
    return palette_data[n];   // fixme: mozna zpet k puvodnimu
  }
  else if (pal_type == 1)
    color[0] = color[1] = color[2] = (float)x;
  else if (pal_type == 2)
    color[0] = color[1] = color[2] = (float)(1.0 - x);
  else
    color[0] = color[1] = color[2] = 1.0;
  return color;

  /*if (n == num_pal_entries) return palette_data[n];
  float s = x - (double) n;
  float t = 1.0 - s;
  static float color[3];
  color[0] = palette_data[n][0] * t + palette_data[n+1][0] * s;
  color[1] = palette_data[n][1] * t + palette_data[n+1][1] * s;
  color[2] = palette_data[n][2] * t + palette_data[n+1][2] * s;
  return color;*/
}


void View::set_num_palette_steps(int num)
{
  if (num < 2) num = 2;
  if (num > 256) num = 256;
  pal_steps = num;
  update_tex_adjust();

  view_sync.enter();
  if (window_id >= 0)
    create_gl_palette();
  view_sync.leave();

  safe_post_redisplay();
}


void View::create_gl_palette()
{
  int i;
  unsigned char palette[256][3];
  for (i = 0; i < pal_steps; i++)
  {
    const float* color = get_palette_color((double) i / pal_steps);
    palette[i][0] = (unsigned char) (color[0] * 255);
    palette[i][1] = (unsigned char) (color[1] * 255);
    palette[i][2] = (unsigned char) (color[2] * 255);
  }
  for (i = pal_steps; i < 256; i++)
    memcpy(palette[i], palette[pal_steps-1], 3);

  if (gl_pallete_tex_id == 0)
    glGenTextures(1, &gl_pallete_tex_id);
  glBindTexture(GL_TEXTURE_1D, gl_pallete_tex_id);
  glTexImage1D(GL_TEXTURE_1D, 0, 3, 256, 0, GL_RGB, GL_UNSIGNED_BYTE, palette);
#ifndef GL_CLAMP_TO_EDGE // fixme: this is needed on Windows
  #define GL_CLAMP_TO_EDGE 0x812F
#endif
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
}


void View::set_palette_filter(bool linear)
{
  view_sync.enter(); //lock to prevent simultaneuous rendering

  pal_filter = linear ? GL_LINEAR : GL_NEAREST;
  
  if (gl_pallete_tex_id == 0)
    glGenTextures(1, &gl_pallete_tex_id);
  glBindTexture(GL_TEXTURE_1D, gl_pallete_tex_id);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, pal_filter);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, pal_filter);
  update_tex_adjust();

  view_sync.leave(); //unlock

  post_redisplay();
}


void View::set_palette(int type)
{
  if (type < 0 || type > 3) error("type can only be 0, 1, 2 or 3.");

  view_sync.enter();
  pal_type = type;
  if (window_id >= 0)
    create_gl_palette();
  view_sync.leave();

  //redisplay
  safe_post_redisplay();
}


void View::update_tex_adjust()
{
  if (pal_filter == GL_LINEAR)
  {
    tex_scale = (double) (pal_steps-1) / 256.0;
    tex_shift = 0.5 / 256.0;
  }
  else
  {
    tex_scale = (double) pal_steps / 256.0;
    tex_shift = 0.0;
  }
}


void View::set_min_max_range(double min, double max)
{
  range_min = min;
  range_max = max;
  range_auto = false;
  if (window_id >= 0) { update_layout(); post_redisplay(); }
}

void View::auto_min_max_range()
{
  range_auto = true;
  if (window_id >= 0) { update_layout(); post_redisplay(); }
}

void View::get_min_max_range(double& min, double& max)
{
  min = range_min;
  max = range_max;
}


void View::draw_text(double x, double y, const char* text, int align)
{
  void* font = GLUT_BITMAP_9_BY_15;
  if (align > -1)
  {
    int width = glutBitmapLength(font, (const unsigned char*) text);
    if (align == 1) x -= width; // align right
               else x -= (double) width / 2; // center
  }
  y += 5; //(double) glutBitmapHeight(font) / 2 - 1;

  glDisable(GL_TEXTURE_1D);
  glDisable(GL_LIGHTING);

  glRasterPos2d((int) (x+0.5), (int) (y+0.5));
  glutBitmapString(font, (const unsigned char*) text);
}


int View::get_text_width(const char* text)
{
  void* font = GLUT_BITMAP_9_BY_15;
  return glutBitmapLength(font, (const unsigned char*) text);
}


void View::draw_help()
{
  set_ortho_projection(true);
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);
  glDisable(GL_TEXTURE_1D);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  const char* text = get_help_text();

  int n = 1;
  for (const char* p = text; *p; p++)
    if (*p == '\n') n++;

  int width = get_text_width(text);
  int height = n * glutBitmapHeight(GLUT_BITMAP_9_BY_15);
  int x = 10, y = 10, b = 6;

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 
  glColor4f(1.0f, 1.0f, 1.0f, 0.65f);
  glBegin(GL_QUADS);
    glVertex2d(x, y+height+2*b);
    glVertex2d(x+width+2*b, y+height+2*b);
    glVertex2d(x+width+2*b, y);
    glVertex2d(x, y);
  glEnd();

  glDisable(GL_BLEND);
  glColor3f(0, 0, 0);
  draw_text(x+b, y+b+7, text);
}


//// screenshot stuff //////////////////////////////////////////////////////////////////////////////

char* View::get_screenshot_file_name()
{
  static char file_name[1024] = { 0 };
  bool got_file_name = false;
  do
  {
    sprintf(file_name, "screen%03d.bmp", screenshot_no);
	  FILE *f = fopen(file_name, "r");
	  if (f == NULL)
	    got_file_name = true;
	  else
  	  fclose(f);
	  screenshot_no++;
  }
  while (!got_file_name);
  return file_name;
}


typedef unsigned int dword;
typedef unsigned short word;

const word BITMAP_ID = 0x4D42;

#pragma pack(1)

struct BitmapFileHeader
{
  word  type;
  dword size;
  word  reserved1;
  word  reserved2;
  dword off_bits;
};

struct BitmapInfoHeader
{
  dword size;
  dword width;
  dword height;
  word  planes;
  word  bit_count;
  dword compression;
  dword size_image;
  dword xdpi;
  dword ydpi;
  dword clr_used;
  dword clr_important;
};


void View::save_screenshot_internal(const char *file_name)
{
  BitmapFileHeader file_header;
  BitmapInfoHeader info_header;

  // alloc memory for pixel data (4 bytes per pixel)
  char* pixels = NULL;
  if ((pixels = (char*) malloc(4 * window_width * window_height)) == NULL)
    error("Could not allocate memory for pixel data");

  // get pixels from framebuffer
#ifdef GL_BGRA_EXT
  glReadPixels(0, 0, window_width, window_height, GL_BGRA_EXT, GL_UNSIGNED_BYTE, pixels);
#else
  glReadPixels(0, 0, window_width, window_height, GL_RGBA, GL_UNSIGNED_BYTE, pixels); // FIXME!!!
  warn("dont have GL_BGRA_EXT format");
#endif

  // opening file for binary writing
  FILE* file = fopen(file_name, "wb");
  if (file == NULL)
    error("Could not open '%s' for writing", file_name);

  // fill in bitmap header
  file_header.type = BITMAP_ID;
  file_header.size = sizeof(BitmapFileHeader) +  sizeof(BitmapInfoHeader) +
                     4 * window_width * window_height;
  file_header.reserved1 = file_header.reserved2 = 0;
  file_header.off_bits = 14 + 40; // length of both headers

  if (fwrite(&file_header, sizeof(file_header), 1, file) != 1)
    error("Error writing bitmap header");

  // fill in bitmap info header
  info_header.size = sizeof(BitmapInfoHeader);
  info_header.width = window_width;
  info_header.height = window_height;
  info_header.planes = 1;
  info_header.bit_count = 32; // 4 bytes per pixel = 32 bits
  info_header.compression = 0;
  info_header.size_image = window_width * window_height * 4;
  info_header.xdpi = 2835; // 72 dpi
  info_header.ydpi = 2835; // 72 dpi
  info_header.clr_used = 0;
  info_header.clr_important = 0;

  if (fwrite(&info_header, sizeof(info_header), 1, file) != 1)
    error("Error writing bitmap header\n");

  // write image pixels
  if (fwrite((GLubyte*) pixels, 1, info_header.size_image, file) != info_header.size_image)
    error("Error writing pixel data\n");

  fclose(file);
  free((void*) pixels);
  info("Saved %s", file_name);
}


void View::save_screenshot(const char* bmpname, bool high_quality)
{
  view_sync.enter();
  if (window_id >= 0) { //set variable neccessary to create a screenshot
    hq_frame = high_quality;
    want_screenshot = true;
    screenshot_filename = bmpname;
  }
  view_sync.leave();

  //request redraw
  safe_post_redisplay(); 
}


void View::save_numbered_screenshot(const char* format, int number, bool high_quality)
{
  char buffer[1000];
  sprintf(buffer, format, number);
  save_screenshot(buffer, high_quality);
}


//// scale /////////////////////////////////////////////////////////////////////////////////////////

int View::measure_scale_labels()
{
  int result = 0;
  for (int i = 0; i <= scale_numticks+1; i++)
  {
    double value = range_min + (double) i * (range_max - range_min) / (scale_numticks+1);
    if (fabs(value) < 1e-8) value = 0.0;
    char text[50];
    sprintf(text, scale_fmt, value);
    int w = get_text_width(text);
    if (w > result) result = w;
  }
  return result;

/*  char text[100];
  sprintf(text, scale_fmt, -0.000123456789123456789);
  return get_text_width(text);*/
}


void View::draw_continuous_scale(char* title, bool righttext)
{
  int i;
  double y0 = scale_y + scale_height;

  set_ortho_projection(true);
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);
  glDisable(GL_TEXTURE_1D);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  // background
  const int b = 5;
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 
  glColor4f(1.0f, 1.0f, 1.0f, 0.65f);
  int rt = righttext ? 0 : labels_width + 8;
  glBegin(GL_QUADS);
  glVertex2d(scale_x - b - rt, y0 + 5 + b);
  glVertex2d(scale_x + scale_width + 8 + labels_width + b - rt, y0 + 5 + b);
  glVertex2d(scale_x + scale_width + 8 + labels_width + b - rt, scale_y - 5 - b);
  glVertex2d(scale_x - b - rt, scale_y - 5 - b);
  glEnd();

  // palette
  glDisable(GL_BLEND);
  glColor3f(0.0f, 0.0f, 0.0f);
  glBegin(GL_QUADS);
  glVertex2d(scale_x, scale_y);
  glVertex2d(scale_x, scale_y + scale_height + 1);
  glVertex2d(scale_x + scale_width + 1, scale_y + scale_height + 1);
  glVertex2d(scale_x + scale_width + 1, scale_y);
  glEnd();

  glEnable(GL_TEXTURE_1D);
  glBindTexture(GL_TEXTURE_1D, gl_pallete_tex_id);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
  glBegin(GL_QUADS);
  glTexCoord1d(tex_scale + tex_shift);
  glVertex2d(scale_x + 1, scale_y + 1);
  glVertex2d(scale_x + scale_width, scale_y + 1);
  glTexCoord1d(tex_shift);
  glVertex2d(scale_x + scale_width, scale_y + scale_height);
  glVertex2d(scale_x + 1, scale_y + scale_height);
  glEnd();

  // focus
  glDisable(GL_TEXTURE_1D);
  if (scale_focused)
  {
    glEnable(GL_BLEND);
    glColor4f(1.0f, 1.0f, 1.0f, 0.3f);
    glBegin(GL_QUADS);
    glVertex2d(scale_x + 1, scale_y + 1);
    glVertex2d(scale_x + scale_width, scale_y + 1);
    glVertex2d(scale_x + scale_width, scale_y + scale_height);
    glVertex2d(scale_x + 1, scale_y + scale_height);
    glEnd();
  }

  // ticks
  glColor3f(0, 0, 0);
  glDisable(GL_BLEND);
  glDisable(GL_LINE_STIPPLE);
  glLineWidth(1.0);
  glBegin(GL_LINES);
  for (i = 0; i < scale_numticks; i++)
  {
    y0 = scale_y + scale_height - (double) (i+1) * scale_height / (scale_numticks+1);
    glVertex2d(scale_x, y0);
    glVertex2d(scale_x + 0.2 * scale_width + 1, y0);
    glVertex2d(scale_x + 0.8 * scale_width, y0);
    glVertex2d(scale_x + scale_width, y0);
  }
  glEnd();

  // labels
  for (i = 0; i <= scale_numticks+1; i++)
  {
    double value = range_min + (double) i * (range_max - range_min) / (scale_numticks+1);
    if (fabs(value) < 1e-8) value = 0.0;
    char text[50];
    sprintf(text, scale_fmt, value);
    y0 = scale_y + scale_height - (double) i * scale_height / (scale_numticks+1);
    if (righttext)
      draw_text(scale_x + scale_width + 8, y0, text);
    else
      draw_text(scale_x - 8, y0, text, 1);
  }

  //if (title != NULL) draw_text(x, y-18, title);
  //if (title != NULL) draw_text(x, y+height+25, title);
}


void View::draw_discrete_scale(int numboxes, const char* boxnames[], const float boxcolors[][3])
{
  set_ortho_projection(true);
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);
  glDisable(GL_TEXTURE_1D);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  // background
  const int b = 5;
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 
  glColor4f(1.0f, 1.0f, 1.0f, 0.65f);
  glBegin(GL_QUADS);
  glVertex2d(scale_x - b, scale_y - b);
  glVertex2d(scale_x - b, scale_y + scale_height + b+1);
  glVertex2d(scale_x + scale_width + b+1, scale_y + scale_height + b+1);
  glVertex2d(scale_x + scale_width + b+1, scale_y - b);
  glEnd();

  // boxes
  glDisable(GL_BLEND);
  int y = scale_y;
  for (int i = 0; i < numboxes; i++)
  {
    glColor3f(0.0, 0.0, 0.0);
    glBegin(GL_QUADS);
    glVertex2d(scale_x, y);
    glVertex2d(scale_x, y + scale_box_height + 1);
    glVertex2d(scale_x + scale_width + 1, y + scale_box_height + 1);
    glVertex2d(scale_x + scale_width + 1, y);
    glEnd();

    const float* color = boxcolors[numboxes-1-i];
    float bcolor[3] = { color[0], color[1], color[2] };
    if (scale_focused) {
      bcolor[0] = color[0]*0.7f + 1.0f*0.3f;
      bcolor[1] = color[1]*0.7f + 1.0f*0.3f;
      bcolor[2] = color[2]*0.7f + 1.0f*0.3f;
    }

    glColor3f(bcolor[0], bcolor[1], bcolor[2]);
    glBegin(GL_QUADS);
    glVertex2d(scale_x+1, y+1);
    glVertex2d(scale_x+1, y + scale_box_height);
    glVertex2d(scale_x + scale_width, y + scale_box_height);
    glVertex2d(scale_x + scale_width, y+1);
    glEnd();

    if ((color[0] + color[1] + color[2]) / 3 > 0.5)
      glColor3f(0, 0, 0);
    else
      glColor3f(1, 1, 1);

    int a = scale_x + scale_width/2;
    int b = y + scale_box_height/2;
    draw_text(a, b, boxnames[numboxes-1-i], 0);
    draw_text(a+1, b, boxnames[numboxes-1-i], 0);

    y += scale_box_height + scale_box_skip;
  }
}


void View::scale_dispatch()
{
  draw_continuous_scale(NULL, !pos_horz);
}


void View::update_layout()
{
  lspace = rspace = labels_width = 0;
  if (b_scale)
  {
    labels_width = scale_fixed_width;
    if (labels_width < 0) labels_width = measure_scale_labels();
    int space = scale_width + 8 + labels_width + margin;
    if (pos_horz == 0)
      { lspace = space;  scale_x = margin; }
    else
      { rspace = space;  scale_x = window_width - margin - scale_width; }

    if (pos_vert == 0)
      scale_y = window_height - margin - scale_height;
    else
      scale_y = margin;
  }

  center_x = ((double) window_width - 2*margin - lspace - rspace) / 2 + margin + lspace;
  center_y = (double) window_height / 2;
}


void View::show_scale(bool show)
{
  b_scale = show;
  if (window_id >= 0) { update_layout(); post_redisplay(); }
}

void View::set_scale_position(int horz, int vert)
{
  pos_horz = horz;
  pos_vert = vert;
  if (window_id >= 0) { update_layout(); post_redisplay(); }
}

void View::set_scale_size(int width, int height, int numticks)
{
  scale_width = width;
  scale_height = height;
  scale_numticks = numticks;
  if (window_id >= 0) { update_layout(); post_redisplay(); }
}

void View::set_scale_format(const char* fmt)
{
  strncpy(scale_fmt, fmt, 19);
  if (window_id >= 0) { update_layout(); post_redisplay(); }
}

void View::fix_scale_width(int width)
{
  scale_fixed_width = width;
  if (window_id >= 0) { update_layout(); post_redisplay(); }
}


#endif // NOGLUT
