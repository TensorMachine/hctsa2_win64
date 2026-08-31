/* tisean_reset.h -- per-program global reset for the TISEAN MEX gateways.
 *
 * A TISEAN program is written to run once and exit. Several of them mutate
 * their own file-scope globals and never restore them: d2 rescales EPSMAX by
 * the data interval and narrows it again during the run, so a second call in
 * the same process starts from the previous call's scaling range and returns
 * a different -- and entirely plausible-looking -- answer.
 *
 * Rather than transcribe each initial value by hand (get one wrong and the
 * failure is silent), the first call snapshots the globals as the compiler
 * left them and every later call restores from that snapshot. The values come
 * from the source, not from this file.
 *
 * Usage in a gateway:
 *
 *     #define D2_GLOBALS(X)   X(double*, series)  X(double, EPSMAX)  ...
 *     TISEAN_DEFINE_RESET(D2_GLOBALS)
 *
 * List every scalar and pointer the program declares at file scope. Large
 * fixed arrays (d2's box[256][256], lyap_r's box[256][256]) are deliberately
 * omitted: each program refills them at the start of every run, and clearing
 * half a megabyte per call would cost more than it buys. Anything omitted must
 * be shown to be write-before-read -- the byte-identity test in
 * HCTSA2_TestTiseanMex is what establishes that.
 */
#ifndef TISEAN_RESET_H
#define TISEAN_RESET_H

#define TS_FIELD(type, name) type name;
#define TS_SNAP(type, name)  ts_init.name = name;
#define TS_REST(type, name)  name = ts_init.name;

#define TISEAN_DEFINE_RESET(LIST)                                    \
    static struct { LIST(TS_FIELD) } ts_init;                        \
    static int ts_have_init;                                         \
    static void reset(void)                                          \
    {                                                                \
        if (!ts_have_init) { LIST(TS_SNAP) ts_have_init = 1; }       \
        else { LIST(TS_REST) }                                       \
    }

#endif /* TISEAN_RESET_H */
