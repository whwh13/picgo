#ifndef MEDYAN_Visual_Common_Hpp
#define MEDYAN_Visual_Common_Hpp

#include "Util/Environment.hpp"

#ifdef PLATFORM_WINDOWS
    #define APIENTRY __stdcall
#endif
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>

namespace medyan::visual {

#ifdef VISUAL
constexpr bool enabled = true;
#else
constexpr bool enabled = false;
#endif

} // namespace medyan::visual

#endif
