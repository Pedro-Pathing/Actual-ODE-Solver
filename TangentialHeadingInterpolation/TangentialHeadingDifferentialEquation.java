package TangentialHeadingInterpolation;

import MathUtil.CubicBezierCurve;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.geometry.euclidean.twod.Vector2D;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.util.FastMath;

import java.awt.geom.Point2D;

public class TangentialHeadingDifferentialEquation implements FirstOrderDifferentialEquations {
    private final double theta;
    private final double v_max; // Maximum velocity (inches per second)
    private final double mass; // Mass of the robot (kg)
    private final double mu_k; // Coefficient of kinetic friction
    private final double c_1;  // Coefficient for lateral penalty
    private final double c_2;  // Coefficient for angular velocity penalty
    private final CubicBezierCurve path;
    private final double angularSpeed;
    private final double feedforward;
    private final double proportional;
    private final double derivative;

    public TangentialHeadingDifferentialEquation(double theta, double vMax, double mass, double muK, double c1, double c2, CubicBezierCurve path, double angularSpeed, double feedForward, double proportional, double derivative) {
        this.theta = theta;
        this.v_max = vMax;
        this.mass = mass;
        this.mu_k = muK;
        this.c_1 = c1;
        this.c_2 = c2;
        this.path = path;
        this.angularSpeed = angularSpeed;
        this.feedforward = feedForward;
        this.proportional = proportional;
        this.derivative = derivative;
    }

    @Override
    public int getDimension() {
        return 2;
    }

    @Override
    public void computeDerivatives(double v, double[] state, double[] derivatives) throws MaxCountExceededException, DimensionMismatchException {
        double t = path.pathInversion(new Point2D.Double(state[0], state[1]));
        Point2D.Double pathPrimes = path.computePathPrimes(t);
        double scalar = v_max - (v_max + c_2 * angularSpeed) * clamp(computeHeadingPower()) - mu_k * mass;
        double heading = FastMath.atan2(pathPrimes.getY(), pathPrimes.getX());

        // Update derivatives
        derivatives[0] = scalar * FastMath.cos(heading); // dx/dt
        derivatives[1] = scalar * FastMath.sin(heading); // dy/dt
    }

    private double computeHeadingPower() {

    }

    private double clamp(double x) {
        return FastMath.max(FastMath.min(x, 1),0);
    }
}