/* Copyright (c) 2014, 2015 Qualcomm Technologies Inc

All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted (subject to the limitations in the disclaimer below) provided that
the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list
of conditions and the following disclaimer.

Redistributions in binary form must reproduce the above copyright notice, this
list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

Neither the name of Qualcomm Technologies Inc nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY THIS
LICENSE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. */

package MathUtil;


import java.util.concurrent.TimeUnit;

/**
 * The {@link ElapsedTime} class provides a simple handy timer to measure elapsed time intervals.
 * The timer does not provide events or callbacks, as some other timers do. Rather, at an application-
 * determined juncture, one can {@link #reset()} the timer. Thereafter, one can query the interval
 * of wall-clock time that has subsequently elapsed by calling the {@link #time()}, {@link #seconds()},
 * or {@link #milliseconds()} methods. The timer has nanosecond internal accuracy. The precision
 * reported by the {@link #time()} method is either seconds or milliseconds, depending on how the
 * timer is initially constructed.
 *
 * This class is thread-safe.
 */

public class ElapsedTime {

    //------------------------------------------------------------------------------------------------
    // Types and constants
    //------------------------------------------------------------------------------------------------

    /**
     * An indicator of the resolution of a timer.
     * @see ElapsedTime#ElapsedTime(Resolution)
     */
    public enum Resolution {
        SECONDS,
        MILLISECONDS
    }

    /** the number of nanoseconds in a second */
    public static final long SECOND_IN_NANO = 1000000000;

    /** the number of nanoseconds in a millisecond */
    public static final long MILLIS_IN_NANO = 1000000;

    //------------------------------------------------------------------------------------------------
    // State
    //------------------------------------------------------------------------------------------------

    protected volatile long nsStartTime;
    protected final double resolution;

    //------------------------------------------------------------------------------------------------
    // Construction
    //------------------------------------------------------------------------------------------------

    public ElapsedTime() {
        reset();
        this.resolution = SECOND_IN_NANO;
    }

    public ElapsedTime(long startTime) {
        this.nsStartTime = startTime;
        this.resolution = SECOND_IN_NANO;
    }

    /**
     * Creates a timer with a resolution of seconds or milliseconds. The resolution
     * affects the units in which the {@link #time()} method reports. The timer is initialized
     * with the current time.
     * @param resolution the resolution of the new timer
     * @see #ElapsedTime()
     */
    public ElapsedTime(Resolution resolution) {
        reset();
        switch (resolution) {
            case SECONDS:
            default:
                this.resolution = SECOND_IN_NANO;
                break;
            case MILLISECONDS:
                this.resolution = MILLIS_IN_NANO;
                break;
        }
    }

    //------------------------------------------------------------------------------------------------
    // Operations
    //------------------------------------------------------------------------------------------------

    protected long nsNow() {
        return System.nanoTime();
    }

    /**
     * Returns the current time on the clock used by the timer
     * @param unit  the time unit in which the current time should be returned
     * @return      the current time on the clock used by the timer
     */
    public long now(TimeUnit unit) {
        return unit.convert(nsNow(), TimeUnit.NANOSECONDS);
    }

    /**
     * Resets the internal state of the timer to reflect the current time. Instantaneously following
     * this reset, {@link #time()} will report as zero.
     * @see #time()
     */
    public void reset() {
        nsStartTime = nsNow();
    }

    /**
     * Returns, in resolution-dependent units, the time at which this timer was last reset.
     * @return the reset time of the timer
     */
    public double startTime() {
        return nsStartTime / resolution;
    }

    /**
     * Returns the time at which the timer was last reset, in units of nanoseconds
     * @return the time at which the timer was last reset, in units of nanoseconds
     */
    public long startTimeNanoseconds() {
        return this.nsStartTime;
    }

    /**
     * Returns the duration that has elapsed since the last reset of this timer.
     * Units used are either seconds or milliseconds, depending on the resolution with
     * which the timer was instantiated.
     * @return time duration since last timer reset
     * @see #ElapsedTime()
     * @see #ElapsedTime(Resolution)
     * @see #seconds()
     * @see #milliseconds()
     */
    public double time() {
        return (nsNow() - nsStartTime) / resolution;
    }

    /**
     * Returns the duration that has elapsed since the last reset of this timer
     * as an integer in the units requested.
     * @param unit  the units in which to return the answer
     * @return time duration since last timer reset
     */
    public long time(TimeUnit unit) {
        return unit.convert(nanoseconds(), TimeUnit.NANOSECONDS);
    }

    /**
     * Returns the duration that has elapsed since the last reset of this timer in seconds
     * @return time duration since last timer reset
     * @see #time()
     */
    public double seconds() {
        return nanoseconds() / ((double)(SECOND_IN_NANO));
    }

    /**
     * Returns the duration that has elapsed since the last reset of this timer in milliseconds
     * @return time duration since last timer reset
     * @see #time()
     */
    public double milliseconds() {
        return seconds() * 1000;
    }

    /**
     * Returns the duration that has elapsed since the last reset of this timer in nanoseconds
     * @return time duration since last timer reset
     * @see #time()
     */
    public long nanoseconds() {
        return (nsNow() - nsStartTime);
    }

    /**
     * Returns the resolution with which the timer was instantiated.
     * @return the resolution of the timer
     */
    public Resolution getResolution() {
        if (this.resolution == MILLIS_IN_NANO)
            return Resolution.MILLISECONDS;
        else
            return Resolution.SECONDS;
    }

    //------------------------------------------------------------------------------------------------
    // Utility
    //------------------------------------------------------------------------------------------------

    private String resolutionStr() {
        if (resolution == SECOND_IN_NANO) {
            return "seconds";
        } else if (resolution == MILLIS_IN_NANO) {
            return "milliseconds";
        } else {
            return "unknown units";
        }
    }

    /**
     * Returns a string indicating the current elapsed time of the timer.
     */
    @Override
    public String toString() {
        return String.format("%1.4f %s", time(), resolutionStr());
    }
}
