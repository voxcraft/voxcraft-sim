#if !defined(UTILS_DEBUG_H)
#define UTILS_DEBUG_H

#define DEBUG_HOST_ENABLED true
#define debugHost(cmd) { if (DEBUG_HOST_ENABLED) {printf(("\n[debugHost] %s(%d): <%s> \033[0;33m"), __FILE__, __LINE__,__FUNCTION__); cmd; printf("\033[0m"); } }
#ifdef DEBUG_HOST_ENABLED
#define debugHostx( var_name, cmd ) { printf(("\n[debugHost] %s(%d): <%s> (%s) \033[0;33m"), __FILE__, __LINE__,__FUNCTION__, var_name); cmd; printf("\033[0m");  }
#else
#define debugHostx( a, b ) ;
#endif

#ifdef __CUDACC__
#define DEBUG_DEV_ENABLED true
#define debugDev(cmd) { if (DEBUG_DEV_ENABLED) {printf(("\n[debugDev] %s(%d): <%s> \033[0;32m"), __FILE__, __LINE__,__FUNCTION__); cmd; printf("\033[0m"); } }
#ifdef DEBUG_DEV_ENABLED
#define debugDevice( var_name, cmd ) { printf(("\n[debugDev] %s(%d): <%s> (%s) \033[0;32m"), __FILE__, __LINE__,__FUNCTION__, var_name); cmd; printf("\033[0m");  }
#else
#define debugDevice( a, b ) ;
#endif
#endif

#endif // UTILS_DEBUG_H
