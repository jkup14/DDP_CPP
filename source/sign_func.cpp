#ifndef SIGN_FUNC_CPP
#define SIGN_FUNC_CPP

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

#endif