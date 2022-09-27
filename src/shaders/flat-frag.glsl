#version 300 es
precision highp float;

uniform vec3 u_Eye, u_Ref, u_Up;
uniform vec2 u_Dimensions;
uniform float u_Time;

const float STEPSIZE = 0.01;
const int MAXSTEPS = 1000;

in vec2 fs_Pos;
in vec4 fs_LightVec;  
out vec4 out_Col;

const vec3 lightPos = vec3(5, 5, 3);

float sdSphere( vec3 p, float s )
{
  return length(p)-s;
}

vec3 lambertShade( vec3 pos, vec3 normal, vec3 diffuse )
{
  // Calculate the diffuse term for Lambert shading
  float diffuseTerm = dot(normalize(normal), normalize(lightPos - pos));
  // Avoid negative lighting values
  // diffuseTerm = clamp(diffuseTerm, 0, 1);

  float ambientTerm = 0.2;

  float lightIntensity = diffuseTerm + ambientTerm;   //Add a small float value to the color multiplier
                                                      //to simulate ambient lighting. This ensures that faces that are not
                                                      //lit by our point light are not completely black.
  // Compute final shaded color
  return diffuse * lightIntensity;
}


void main() {

  vec3 eye = u_Eye;
  vec3 forward = u_Ref - u_Eye;
  vec3 up = u_Up;
  vec3 right = normalize(cross(forward,up));

  float f = 0.05 * distance(u_Eye, u_Ref);
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

  do {
    d = sdSphere(pos, 1.0);
    if (d < 0.01)
    {
      col = lambertShade(pos, pos, vec3(1.0));
      break;
    }
    pos += d*dir;
    ++itr;
  } while(itr < MAXSTEPS);

  // vec4 surfaceColor = vec4(0.0);
  // vec4 boxColor = vec4(0.0);
  // vec4 background = vec4(0.9);
  // vec4 wallColor = background;

  // vec3 p = eye;
  out_Col = vec4(col, 1.0);
  // out_Col = vec4(0.5 * (fs_Pos + vec2(1.0)), 0.5 * (sin(u_Time * 3.14159 * 0.01) + 1.0), 1.0);
}
