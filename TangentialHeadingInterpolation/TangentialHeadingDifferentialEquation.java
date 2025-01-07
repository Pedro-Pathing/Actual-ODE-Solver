package TangentialHeadingInterpolation;

import MathUtil.CubicBezierCurve;
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

    public TangentialHeadingDifferentialEquation(double theta, double vMax, double mass, double muK, double c1, double c2, CubicBezierCurve path) {
        this.theta = theta;
        this.v_max = vMax;
        this.mass = mass;
        this.mu_k = muK;
        this.c_1 = c1;
        this.c_2 = c2;
        this.path = path;
    }

    @Override
    public int getDimension() {
        return 2;
    }

    @Override
    public void computeDerivatives(double t, double[] state, double[] derivatives) throws MaxCountExceededException, DimensionMismatchException {
        Point2D.Double position = new Point2D.Double(state[0], state[1]);
        double pathParam = path.pathInversion(position);
        Point2D.Double pathPrime = path.computePathPrimes(pathParam);
        Vector2D pathPrimes = new Vector2D(pathPrime.getX(), pathPrime.getY());

        // Compute heading (H)
        double heading = computeHeading(pathParam);

        // Compute friction term
        double friction = computeFriction(heading);

        // Compute scalar (tangential component of velocity)
        double scalar = computeScalar(pathPrimes, heading, friction);

        // Update derivatives
        derivatives[0] = scalar * pathPrimes.getX(); // dx/dt
        derivatives[1] = scalar * pathPrimes.getY(); // dy/dt
    }

    // Method for computing heading (H(t))
    private double computeHeading(double pathParam) {
        Point2D.Double pathPrime = path.computePathPrimes(pathParam);
        return FastMath.atan2(pathPrime.getY(), pathPrime.getX());
    }

    // Method for computing friction (Fk(theta))
    private double computeFriction(double heading) {
        double angularVelocity = computeAngularVelocity(); // Placeholder for actual angular velocity computation
        return mu_k * mass + c_1 * FastMath.abs(FastMath.sin(heading)) + c_2 * FastMath.pow(angularVelocity, 2);
    }

    // Method for computing scalar (tangential velocity adjustment)
    private double computeScalar(Vector2D pathPrimes, double heading, double friction) {
        double dotProduct = pathPrimes.dotProduct(new Vector2D(FastMath.cos(theta), FastMath.sin(theta)));
        double normPathPrimes = pathPrimes.getNorm();
        return (v_max - friction) * (dotProduct / normPathPrimes);
    }

    // Placeholder method for angular velocity (dtheta/dt)
    private double computeAngularVelocity() {
        // Implement angular velocity computation if needed
        return 0; // Replace with actual value
    }
}