#ifndef MEDYAN_Visual_Shader_Hpp
#define MEDYAN_Visual_Shader_Hpp

#include "Util/Io/Log.hpp"
#include "Visual/Common.hpp"

namespace medyan::visual {

// Managing the shader compilation process
//
// Note:
//   - Creation/modification/destruction of objects of this type must be within
//     an opengl context
class Shader {
public:

    Shader(const char* vertexShaderSrc, const char* fragmentShaderSrc) {
        init_(vertexShaderSrc, fragmentShaderSrc);
    }

    ~Shader() {
        destroy_();
    }

    auto id() const { return id_; }

    void reset(const char* vertexShaderSrc, const char* fragmentShaderSrc) {
        destroy_();
        init_(vertexShaderSrc, fragmentShaderSrc);
    }

    void setBool(const char* name, bool val) const {
        glUniform1i(glGetUniformLocation(id_, name), val);
    }
    void setFloat(const char* name, float f) const {
        glUniform1f(glGetUniformLocation(id_, name), f);
    }
    void setMat3(const char* name, const glm::mat3& mat) const {
        glUniformMatrix3fv(glGetUniformLocation(id_, name), 1, GL_FALSE, glm::value_ptr(mat));
    }
    void setMat4(const char* name, const glm::mat4& mat) const {
        glUniformMatrix4fv(glGetUniformLocation(id_, name), 1, GL_FALSE, glm::value_ptr(mat));
    }
    void setVec3(const char* name, const glm::vec3& vec) const {
        glUniform3fv(glGetUniformLocation(id_, name), 1, glm::value_ptr(vec));
    }

private:
    unsigned int id_;

    void init_(const char* vertexShaderSrc, const char* fragmentShaderSrc) {
        int success;
        char infoLog[512];

        unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vertexShader, 1, &vertexShaderSrc, NULL);
        glCompileShader(vertexShader);
        glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
        if (!success) {
            glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
            log::error("Vertex shader compile failed: {}", infoLog);
        }

        unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fragmentShader, 1, &fragmentShaderSrc, NULL);
        glCompileShader(fragmentShader);
        glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
        if (!success) {
            glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog);
            log::error("Fragment shader compile failed: {}", infoLog);
        }

        id_ = glCreateProgram();
        glAttachShader(id_, vertexShader);
        glAttachShader(id_, fragmentShader);
        glLinkProgram(id_);
        glGetProgramiv(id_, GL_LINK_STATUS, &success);
        if (!success) {
            glGetProgramInfoLog(id_, 512, NULL, infoLog);
            log::error("Shader program link failed: {}", infoLog);
        }

        glDetachShader(id_, vertexShader);
        glDetachShader(id_, fragmentShader);
        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);

    }

    void destroy_() const {
        glDeleteProgram(id_);
    }
};

} // namespace medyan::visual

#endif
