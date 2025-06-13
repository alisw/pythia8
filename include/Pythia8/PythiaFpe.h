// PythiaFpe.h is a part of the PYTHIA event generator.
// Copyright (C) 2025 Christian Bierlich and Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Definition of behaviour for debugging flag -DGCCFPDEBUG. This
// snippet should, when enabled, cause programs to fail at runtime if
// they contain a floating point exception. Note that underflow and
// inexact are not caught, see comment in raisefpe() to enable
// them. This is not recommended, as FE_UNDERFLOW would also catch
// things like small exponents, and FE_INEXACT would catch things like
// 1.0 / 3.0. Expected output is "Caught SIGFPE (Floating Point
// Exception)" followed by a stack trace. The program will exit with
// error. Note that this behaviour is both compiler dependent and
// platform dependent. It works only with gcc on x86/x86_64 (with
// SSE). On other platforms (ARM, PPC, etc.), SSE intrinsics won't be
// valid. This includes x86_64 emulated on ARM e.g. in docker.

// Run ./configure with --obj-common=-DGCCFPDEBUG to enable.
// For better stack trace symbols add -O0 -g -rdynamic:
// --obj-common='-g -O0 -rdynamic -DGCCFPDEBUG'
// Consider also -lexecinfo (if you see "undefined reference to
// backtrace...") and -fno-omit-frame-pointer.
// Should never be combined with compile flags such as:
// -ffast-math -Ofast -fno-trapping-math
// or other fast math flags, as they may optimize away the
// trapping.

#ifndef Pythia8_PythiaFpe_H
#define Pythia8_PythiaFpe_H

// Check the GCCFPDEBUG flag.
#ifdef GCCFPDEBUG

// Catch compilation on ARM platforms.
#ifdef __aarch64__
#error "GCCDEBUG unsupported on ARM64. Disable -DGCCDEBUG or compile on x86."
#endif

#ifndef __ENABLE_FP_DEBUG__
#define __ENABLE_FP_DEBUG__
#include <csignal>     // Provides sigaction, siginfo_t.
#include <fenv.h>      // Provides feenableexcept().
#include <xmmintrin.h> // Provides _MM_GET/SET_EXCEPTION_MASK().
#include <execinfo.h>  // Provides backtrace(), backtrace_symbols_fd().
#include <unistd.h>    // Provides STDERR_FILENO.

// Implement a signal handler. This will catch the action raised by
// the FPE trigger, print a stack trace and exit. This makes it easier
// to find the location in the code where the FPE happened. Handle
// the signal.
static void fpeSignalHandler(int sig, siginfo_t* info, void* context) {

  // Suppress compiler warnings from -Wunused-parameter.
  (void)sig;
  (void)info;
  (void)context;

  // Print an error message. Avoid C++ I/O for async safety.
  fprintf(stderr, "\n*************************************************\n");
  fprintf(stderr,  "** Caught SIGFPE (Floating Point Exception)    **\n");
  fprintf(stderr,  "** Printing stack trace (compile with -O0 -g)  **\n");
  fprintf(stderr,  "** For better symbols, also consider -rdynamic.**\n");
  fprintf(stderr,  "*************************************************\n");

  // Obtain a backtrace.
  void* buffer[32];
  int n = backtrace(buffer, 32);

  // Print the backtrace.
  backtrace_symbols_fd(buffer, n, STDERR_FILENO);
  fprintf(stderr,  "*************************************************\n");

  // Exit with error. Use low-level _exit to avoid calling destructors.
  _exit(1);

}

// Setup a handler to catch the signal allowing an action before the
// program exits.
static void __attribute__((constructor)) setupFpeHandler() {

  // Setup a sigaction for the handler.
  struct sigaction sa;
  // Define the action.
  sa.sa_sigaction = fpeSignalHandler;
  // Empty the handler.
  sigemptyset(&sa.sa_mask);
  // Set the action we want: siginfo_t
  sa.sa_flags = SA_SIGINFO;
  // Install the handler
  sigaction(SIGFPE, &sa, nullptr);

}

// Enable raising of the FPE by unmasking FPU exceptions. This
// includes both x87 FPU (feenableexcept) and SSE exceptions.
static void __attribute__((constructor)) raisefpe() {

  // Enable x87 FPU exceptions. To catch all exceptions, add
  // FE_UNDERFLOW and FE_INEXACT.
  feenableexcept (FE_DIVBYZERO | FE_OVERFLOW | FE_INVALID);

  // Enable SSE exceptions. To catch all exceptions, add
  // _MM_MASK_UNDERFLOW and _MM_MASK_INEXACT.
  unsigned int cw = _MM_GET_EXCEPTION_MASK();
  cw &= ~(_MM_MASK_DIV_ZERO | _MM_MASK_INVALID | _MM_MASK_OVERFLOW);
  _MM_SET_EXCEPTION_MASK(cw);

}

//==========================================================================

#endif // GCCFPDEBUG

#endif // __ENABLE_FP_DEBUG__

#endif // Pythia8_PythiaFpe_H
