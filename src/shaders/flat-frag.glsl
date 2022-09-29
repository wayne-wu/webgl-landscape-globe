#version 300 es
precision highp float;

uniform vec3 u_Eye, u_Ref, u_Up;
uniform vec2 u_Dimensions;
uniform float u_Time;

uniform mat4 u_ViewProj;

const float STEPSIZE = 0.01;
const int MAXSTEPS = 1000;
const float EPS = 0.001;
const float DISCRETIZE_NUM = 4.0;

const float SPHERE_RADIUS = 1.5;
const float TERRAIN_FREQ = 2.0;

const vec3 SKYCOLOR = vec3(0.47, 0.66, 0.82);

const vec3 KEYLIGHT_POS = vec3(15, 15, 10);
const vec3 KEYLIGHT = vec3(1.0, 1.0, 0.9) * 1.5;
const vec3 FILLLIGHT_POS = vec3(0, -5, 0);
const vec3 FILLLIGHT = vec3(1.0, 1.0, 0.9) * 1.0;
const vec3 BACKLIGHT_POS = vec3(0, 5, -5);
const vec3 BACKLIGHT = vec3(1.0, 1.0, 0.9) * 1.0;

in vec2 fs_Pos;
in vec4 fs_LightVec;  
out vec4 out_Col;


float ease_in_quadratic(float t) 
{
    return t * t;
}

float bias(float b, float t) 
{
    return pow(t, log(b) / log(0.5));
}

float gain(float g, float t) 
{
    if (t < 0.5)
        return bias(1.0-g, 2.0*t)/2.0;
    else
        return 1.0 - bias(1.0-g, 2.0 - 2.0*t)/2.0;
}

float impulse(float k, float x)
{
    float h = k*x;
    return h * exp(1.0-h);
}

float map(float value, float min1, float max1, float min2, float max2) {
  return min2 + (value - min1) * (max2 - min2) / (max1 - min1);
}

// Hash functions are taken from IQ's shadertoy examples
float hash (in vec2 st) {
    return fract(sin(dot(st.xy,
                  vec2(12.9898,78.233)))
                 * 43758.5453123);
}

vec3 hash3( vec3 p )
{
	p = vec3( dot(p,vec3(127.1,311.7, 74.7)),
			  dot(p,vec3(269.5,183.3,246.1)),
			  dot(p,vec3(113.5,271.9,124.6)));

	return -1.0 + 2.0*fract(sin(p)*43758.5453123);
}

float trilinear(float a, float b, float c, float d, 
                float e, float f, float g, float h, vec3 u)
{
    return mix(mix(mix(a, b, u.x), mix(c, d, u.x), u.y), 
                mix(mix(e, f, u.x), mix(g, h, u.x), u.y), u.z);
}

float bilinear(float a, float b, float c, float d, vec2 u)
{
  return mix(mix(a, b, u.x), mix(c, d, u.x), u.y);
}

vec2 cubic(vec2 t)
{
    return t*t*(3.0-2.0*t);
}

vec2 quintic(vec2 t)
{
    return t * t * t * (t * (t * 6.0 - 15.0) + 10.0);
}

// float grad(vec3 2, vec3 2, vec2 inc)
// {
//     return dot(hash3(i + inc), f - inc);
// }

float noise( in vec2 x )
{
    vec2 i = floor(x);
    vec2 u = fract(x);
    u = quintic(u);

    float a = hash(i + vec2(0,0));
    float b = hash(i + vec2(1,0));
    float c = hash(i + vec2(0,1));
    float d = hash(i + vec2(1,1));

    return bilinear(a, b, c, d, u);
}

float fbm(in vec2 pos)
{
    float total = 0.f;
    float amplitudeSum = 0.f;

    for (int i = 0; i < 6; i++)
    {
        float frequency = pow(2.0f, float(i));
        float amplitude = pow(0.4f, float(i));
        
        amplitudeSum += amplitude;

        total += amplitude*noise(frequency*pos*1.0);
    }

    return total/amplitudeSum;
}


float sdBox(vec3 p, vec3 b)
{
    vec3 d = abs(p) - b;
    return min(max(d.x, max(d.y, d.z)), 0.0) + length(max(d, 0.0));
}

float sdBBox(vec3 p)
{
    return sdBox(p, vec3(2.0, 2.0, 2.0));
}

float sdSphere( vec3 p, float s )
{
  return length(p)-s;
}

float sdHeight( in vec3 p )
{
  float h = fbm(TERRAIN_FREQ*p.xz);
  h = map(h, 0.0, 1.0, -0.8*SPHERE_RADIUS, 0.5*SPHERE_RADIUS);
  return p.y - h;
}

float discr( float x )
{
  float w = 1.0/DISCRETIZE_NUM;
  return floor(x/w)*w;
}

float shadowP (in vec3 pos, in vec3 dir)
{
  pos += dir * 0.1;
  for (int i = 0; i < 50; ++i)
  {
    float d = sdHeight( pos );
    if (d < EPS)
      return 0.0;
    pos += dir * d;
  }
  return 1.0;
}

vec3 shadeToon( vec3 pos, vec3 normal, vec3 albedo )
{
  vec3 n = normalize(normal);
  // Calculate the diffuse term for Lambert shading
  float d1 = discr(clamp(dot(n, normalize(KEYLIGHT_POS - pos)), 0.0, 1.0));
  float d2 = discr(clamp(dot(n, normalize(FILLLIGHT_POS - pos)), 0.0, 1.0));
  float d3 = discr(clamp(dot(n, normalize(BACKLIGHT_POS - pos)), 0.0, 1.0));

  float ambientTerm = 0.2;

  float lightIntensity = d1 + d2 + d3 + ambientTerm;

  return albedo * lightIntensity;
}

vec3 shadeLambert( vec3 pos, vec3 normal, vec3 albedo )
{
  vec3 n = normalize(normal);
  // Calculate the diffuse term for Lambert shading

  vec3 col = albedo * 
            KEYLIGHT * 
            max(0.0, dot(n, normalize(KEYLIGHT_POS - pos)));
            // shadowP(pos, normalize(KEYLIGHT_POS - pos));
  col += albedo * FILLLIGHT * clamp(dot(n, normalize(FILLLIGHT_POS - pos)), 0.0, 1.0);
  col += albedo * BACKLIGHT * clamp(dot(n, normalize(BACKLIGHT_POS - pos)), 0.0, 1.0);

  // float ambientTerm = 0.2;

  // float lightIntensity = d1 + d2 + d3 + ambientTerm;

  return col;
}

vec3 skyColor( vec2 uv )
{
  float t = map(uv.y, -1.0, 1.0, 0.0, 1.0);
  return mix(vec3(0.0), SKYCOLOR, t);
}

vec3 heightColor( float h )
{
  vec3 color = vec3(0.0);
  if (h > 0.1)
    color = vec3(0.94, 0.95, 0.93);
  else if (h > -0.5)
    color = vec3(0.44, 0.47, 0.27);
  else
    color = vec3(0.46, 0.38, 0.33);
  return color;
}

bool hitBSphere( inout vec3 p, in vec3 dir)
{
  for (int i = 0; i < 50; i++)
  {
    float d = sdSphere(p, SPHERE_RADIUS);
    if (d < EPS)
      return true;
    p += dir * d;
  }
  return false;
}

bool hitTerrain( inout vec3 p, in vec3 dir, inout vec3 color )
{
  for (int i = 0; i < MAXSTEPS; i++)
  {
      // When we hit very close to the surface
      float d = sdHeight(p);
      if (d < EPS)
      {
          // Calculate normal
          vec3 dx = vec3(EPS, 0, 0);
          vec3 dz = vec3(0, 0, EPS);
          vec3 n = vec3(sdHeight(p - dx) - sdHeight(p - dx), 
                        2.0*EPS,
                        sdHeight(p - dz) - sdHeight(p + dz));
          n = normalize(n);

          // if (abs(dot(dir, n))<0.1)
          //   color = vec3(0.0);
          // else
          color += shadeLambert(p, n, heightColor(p.y));
          return true;
      }

      // Checks if the ray hits the bounding box from the inside
      if (sdBBox(p) >= EPS)
      {
          dir = refract(dir, normalize(p), 1.0/1.2);
          p += dir*0.5;
          vec4 uv = u_ViewProj*vec4(p, 1.0);

          color += skyColor(vec2(uv.x, uv.y));
          return false;
      }

      p += dir * STEPSIZE;
  }
  return true;
}

void main() {

  vec3 eye = u_Eye;
  vec3 forward = normalize(u_Ref - u_Eye);
  vec3 up = u_Up;
  vec3 right = normalize(cross(forward,up));

  float f = 0.5 * distance(u_Eye, u_Ref);
  float u = gl_FragCoord.x * 2.0 / u_Dimensions.x - 1.0;
  float v = gl_FragCoord.y * 2.0 / u_Dimensions.y - 1.0;

  float aspectRatio = u_Dimensions.x / u_Dimensions.y;
  right *= aspectRatio;

  // ray's world position
  vec3 pos = eye + right * u + up * v + forward * f;
  // ray's direction
  vec3 dir = normalize(pos - eye);

  int itr = 0;
  float d = 0.0;

  vec3 col = vec3(0.0);
  if (hitBSphere(pos, dir))
  {
    if(sdHeight(pos) < EPS)
      col = shadeLambert(pos, normalize(pos), heightColor(pos.y));
    else
    {
      dir = refract(dir, normalize(pos), 1.0/1.2);
      hitTerrain(pos, dir, col);
    }
  }
  else 
  {
    col = skyColor(vec2(u, v));
  }

  // do {
  //   d = sdSphere(pos, 1.0);
  //   if (d < EPS)
  //   {
  //     col = shadeLambert(pos, pos, vec3(1.0));
  //     break;
  //   }
  //   pos += d*dir;
  //   ++itr;
  // } while(itr < MAXSTEPS);

  // vec4 surfaceColor = vec4(0.0);
  // vec4 boxColor = vec4(0.0);
  // vec4 background = vec4(0.9);
  // vec4 wallColor = background;

  // vec3 p = eye;
  out_Col = vec4(col, 1.0);
  // out_Col = vec4(0.5 * (fs_Pos + vec2(1.0)), 0.5 * (sin(u_Time * 3.14159 * 0.01) + 1.0), 1.0);
}
