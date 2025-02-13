package TangentialHeadingInterpolation;

import MathUtil.CubicBezierCurve;
import MathUtil.ElapsedTime;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.geometry.euclidean.twod.Vector2D;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.util.FastMath;

import java.awt.geom.Point2D;

public class TangentialHeadingDifferentialEquation implements FirstOrderDifferentialEquations {
    private final double v_max; // Maximum velocity (inches per second)
    private final double mass; // Mass of the robot (kg)
    private final double mu_k; // Coefficient of kinetic friction
    private final double c_2;  // Coefficient for angular velocity penalty
    private final CubicBezierCurve path;
    private final double angularVelocity;
    private final double f_heading;
    private final double d_heading;
    private final double p_heading;
    private ElapsedTime timer = new ElapsedTime();
    private double previousHeadingEstimate;
    private double previousHeadingDerivEstimate;

    public TangentialHeadingDifferentialEquation(double vMax, double mass, double muK, double c2, CubicBezierCurve path, double angularVelocity, double fHeading, double dHeading, double pHeading) {
        this.v_max = vMax;
        this.mass = mass;
        this.mu_k = muK;
        this.c_2 = c2;
        this.path = path;
        previousHeadingEstimate = path.theta(0);
        this.angularVelocity = angularVelocity;
        f_heading = fHeading;
        d_heading = dHeading;
        p_heading = pHeading;
        previousHeadingDerivEstimate = 0;
        timer.reset();
    }

    @Override
    public int getDimension() {
        return 2;
    }

    @Override
    public void computeDerivatives(double v, double[] state, double[] derivatives) throws MaxCountExceededException, DimensionMismatchException {
        double currentTime = timer.milliseconds();
        timer.reset();
        double t = path.pathInversion(new Point2D.Double(state[0], state[1]));
        Point2D.Double pathPrime = path.computePathPrimes(t);
        Vector2D pathPrimes = new Vector2D(pathPrime.getX(), pathPrime.getY());
        double headingTarget = FastMath.atan2(pathPrime.getY(), pathPrime.getX());
        double headingDerivativeApprox = (headingTarget - previousHeadingEstimate)/currentTime/1000;
        double headingSecondDerivApprox = (headingDerivativeApprox - previousHeadingDerivEstimate)/currentTime/1000;
        double headingPower = computeHeadingPower(headingDerivativeApprox, headingSecondDerivApprox);
        double scalar = (1 - headingPower) * (v_max - pathPrimes.getNorm() * getFriction(headingPower)) / pathPrimes.getNormSq();

        // Update derivatives
        derivatives[0] = scalar * pathPrimes.getX(); // dx/dt
        derivatives[1] = scalar * pathPrimes.getY(); // dy/dt
        previousHeadingEstimate = headingTarget;
        previousHeadingDerivEstimate = headingDerivativeApprox;
    }

    public double computeHeadingPower(double headingDerivativeApprox, double headingSecondDerivativeApprox) {
        return clamp(f_heading - p_heading * FastMath.abs(headingDerivativeApprox) - d_heading * FastMath.abs(headingSecondDerivativeApprox));
    }

    private double clamp(double x) {
        return FastMath.max(FastMath.min(x, 1),0);
    }

    private double getFriction(double headingPower) {
        return mu_k * mass * 9.81 + c_2*FastMath.abs(angularVelocity * headingPower);
    }
}