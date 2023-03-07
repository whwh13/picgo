#ifndef MEDYAN_Visual_KeyMapping_hpp
#define MEDYAN_Visual_KeyMapping_hpp

#include <array>
#include <string>
#include <unordered_map>
#include <utility> // std::move

#include "utility.h" // std::combine_hash
#include "Visual/Common.hpp"

namespace medyan::visual {

struct KeyComb {
    int glfwKey = GLFW_KEY_UNKNOWN;
    bool ctrl = false;
    bool alt = false;
    bool shift = false;

    bool isGlfwKeyValid() const {
        return (glfwKey >= GLFW_KEY_0 && glfwKey <= GLFW_KEY_9)
            || (glfwKey >= GLFW_KEY_A && glfwKey <= GLFW_KEY_Z)
            || (glfwKey == GLFW_KEY_SPACE);
    }
};

// Default it in C++20.
inline bool operator==(const KeyComb& k1, const KeyComb& k2) {
    return k1.glfwKey == k2.glfwKey && k1.ctrl == k2.ctrl && k1.alt == k2.alt && k1.shift == k2.shift;
}

struct KeyCombHash {
    std::size_t operator()(const KeyComb& k) const {
        std::size_t res = 0;
        std::hash_combine(res, k.glfwKey);
        std::hash_combine(res, k.ctrl);
        std::hash_combine(res, k.alt);
        std::hash_combine(res, k.shift);
        return res;
    }
};

inline std::string display(const KeyComb& kc) {
    if(!kc.isGlfwKeyValid()) {
        return "";
    }
    else {
        const char* keyName;
        switch(kc.glfwKey) {
            case GLFW_KEY_SPACE: keyName = "space"; break;
            default:             keyName = glfwGetKeyName(kc.glfwKey, 0); break;
        }
        std::string res;
        res += kc.ctrl  ? "ctrl+" : "";
        res += kc.alt   ? "alt+" : "";
        res += kc.shift ? "shift+" : "";
        res += keyName  ? keyName : "?";
        return res;
    }
}

enum class KeyCallbackFunction {
    takeSnapshot,
    takeSnapshotAll,
    toggleGui,
    togglePlayPause,
    controlMouseSetRotate,
    controlMouseSetPan,
    last_
};
constexpr auto text(KeyCallbackFunction val) {
    switch(val) {
        case KeyCallbackFunction::takeSnapshot:          return "Take snapshot";
        case KeyCallbackFunction::takeSnapshotAll:       return "Take snapshot for all frames";
        case KeyCallbackFunction::toggleGui:             return "Toggle GUI";
        case KeyCallbackFunction::togglePlayPause:       return "Toggle play/pause";
        case KeyCallbackFunction::controlMouseSetRotate: return "Mouse mode rotate";
        case KeyCallbackFunction::controlMouseSetPan:    return "Mouse mode pan";
        default:                                         return "";
    }
}


struct KeyMapping {
    static KeyMapping createDefault() {
        KeyMapping res;
        res.update(KeyComb { GLFW_KEY_F }, KeyCallbackFunction::takeSnapshot);
        res.update(KeyComb { GLFW_KEY_P }, KeyCallbackFunction::takeSnapshotAll);
        res.update(KeyComb { GLFW_KEY_G }, KeyCallbackFunction::toggleGui);
        res.update(KeyComb { GLFW_KEY_SPACE }, KeyCallbackFunction::togglePlayPause);
        res.update(KeyComb { GLFW_KEY_R }, KeyCallbackFunction::controlMouseSetRotate);
        res.update(KeyComb { GLFW_KEY_T }, KeyCallbackFunction::controlMouseSetPan);
        return res;
    }

    std::unordered_map<KeyComb, KeyCallbackFunction, KeyCombHash> keyCallbackMap;
    std::array< KeyComb, underlying(KeyCallbackFunction::last_) > keyBindingList {};

    // If keyComb is valid:
    //   - clear previous mapping involving this key if any.
    //   - clear previous mapping involving this callback if any.
    //   - add mapping between the key and the callback.
    // Otherwise (if keyComb is invalid):
    //   - clear previous mapping involving this callback if any.
    void update(const KeyComb& keyComb, KeyCallbackFunction func) {
        const auto clearKey = [&keyComb, this] {
            // Find the function it maps to.
            auto it = keyCallbackMap.find(keyComb);
            if(it != keyCallbackMap.end()) {
                // Set to default (invalid) state.
                keyBindingList[underlying(it->second)] = KeyComb {};
                // Clear key in map.
                keyCallbackMap.erase(it);
            }
        };
        const auto clearCallback = [func, this] {
            // Find the keyComb it uses.
            auto& keyComb = keyBindingList[underlying(func)];
            // Clear usage of this key, if any.
            keyCallbackMap.erase(keyComb);
            // Clear keyComb of func.
            keyComb = KeyComb {};
        };
        const auto addMapping = [&, this] {
            keyCallbackMap[keyComb] = func;
            keyBindingList[underlying(func)] = keyComb;
        };

        if(keyComb.isGlfwKeyValid()) {
            // Valid keyComb.
            clearKey();
            clearCallback();
            addMapping();
        }
        else {
            // Invalid keyComb.
            clearCallback();
        }
    }
};



} // namespace medyan::visual

#endif
