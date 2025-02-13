package LinearHeadingInterpolation;

import MathUtil.CubicBezierCurve;
import MathUtil.ElapsedTime;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.genetics.FixedElapsedTime;
import org.apache.commons.math3.geometry.euclidean.twod.Vector2D;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.util.FastMath;

import java.awt.geom.Point2D;
import java.util.Timer;

public class LinearHeadingDifferentialEquation implements FirstOrderDifferentialEquations {
    private final double deltaTheta;
    private final double theta_0;
    private final double v_max; //inches per second
    private final double mass; //kg
    private final double mu_k;
    private final double c_1;
    private final double c_2;
    private final double f_heading;
    private final double d_heading;
    private final double p_heading;
    private final double angularVelocity;
    private final CubicBezierCurve path;
    private double previousHeadingEstimate;
    private double previousHeadingDerivEstimate;
    private ElapsedTime timer = new ElapsedTime();

    public LinearHeadingDifferentialEquation(double deltaTheta, double theta0, double vMax, double mass, double muK, double c1, double c2, double fHeading, double dHeading, double pHeading, double angularVelocity, CubicBezierCurve path) {
        this.deltaTheta = deltaTheta;
        theta_0 = theta0;
        v_max = vMax;
        this.mass = mass;
        mu_k = muK;
        c_1 = c1;
        c_2 = c2;
        f_heading = fHeading;
        d_heading = dHeading;
        p_heading = pHeading;
        this.angularVelocity = angularVelocity;
        this.path = path;
        previousHeadingEstimate = theta0;
        previousHeadingDerivEstimate = 0;
        timer.reset();
    }

    @Override
    public int getDimension() {
        return 2;
    }

    @Override
    public void computeDerivatives(double v, double[] doubles, double[] doubles1) throws MaxCountExceededException, DimensionMismatchException {
        double currentTime = timer.milliseconds();
        timer.reset();
        double t = path.pathInversion(new Point2D.Double(doubles[0], doubles[1]));
        Point2D.Double pathPrime = path.computePathPrimes(t);
        Vector2D pathPrimes = new Vector2D(pathPrime.getX(), pathPrime.getY());
        double theta = deltaTheta * t + theta_0;
        Vector2D angleVector = new Vector2D(FastMath.cos(theta), FastMath.sin(theta));
        double headingTarget = FastMath.atan2(pathPrime.getY(), pathPrime.getX());
        double headingDerivativeApprox = (headingTarget - previousHeadingEstimate)/currentTime/1000;
        double headingSecondDerivApprox = (headingDerivativeApprox - previousHeadingDerivEstimate)/currentTime/1000;
        double headingPower = computeHeadingPower(headingDerivativeApprox, headingSecondDerivApprox);
        double scalar = pathPrimes.dotProduct(angleVector)*v_max - pathPrimes.getNorm()*getFriction(theta, headingTarget, headingPower);
        scalar = scalar * (1 - headingPower) / pathPrimes.getNormSq();

        doubles[0] = scalar * pathPrimes.getX();
        doubles[1] = scalar * pathPrimes.getY();
        previousHeadingEstimate = headingTarget;
        previousHeadingDerivEstimate = headingDerivativeApprox;
    }

    public double computeHeadingPower(double headingDerivativeApprox, double headingSecondDerivativeApprox) {
        return clamp(f_heading - p_heading * FastMath.abs(headingDerivativeApprox) - d_heading * FastMath.abs(headingSecondDerivativeApprox));
    }

    private double clamp(double x) {
        return FastMath.max(FastMath.min(x, 1),0);
    }

    private double getFriction(double theta, double theta_desired, double headingPower) {
        return mu_k * mass * 9.81 + c_1 * FastMath.sin(FastMath.abs(theta - theta_desired))+c_2*FastMath.abs(angularVelocity * headingPower);
    }
}
