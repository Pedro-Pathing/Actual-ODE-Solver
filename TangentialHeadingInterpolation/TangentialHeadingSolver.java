package TangentialHeadingInterpolation;

import MathUtil.CubicBezierCurve;
import MathUtil.SolutionPoints;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.ExpandableStatefulODE;
import org.apache.commons.math3.ode.nonstiff.GillIntegrator;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.OptimizationData;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultiStartMultivariateOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import org.apache.commons.math3.random.RandomVectorGenerator;
import org.apache.commons.math3.util.FastMath;

import java.awt.geom.Point2D;

public class TangentialHeadingSolver {
    // Parameters
    private final double v_max; // Maximum velocity (inches per second)
    private final double mass; // Mass of the robot (kg)
    private final double mu_k; // Coefficient of kinetic friction
    private final double c_1;  // Coefficient for lateral penalty
    private final double c_2;  // Coefficient for angular velocity penalty
    private final Point2D p_0; // Start point
    private final Point2D p_3; // End point
    private final double theta_final; // Final heading angle
    private final boolean isBlueAlliance;
    private final double theta_initial; // Initial heading angle
    private final double angularVelocity; // Maximum angular velocity
    private final double boundaryTolerance; // Distance from boundaries
    private final double submersibleTolerance; // Tolerance for submersibles
    private final double robotWidth; // Width of the robot
    private final double robotHeight; // Height of the robot

    // Integrator and Output Model
    private final GillIntegrator integrator;
    private final ContinuousOutputModel outputModel;

    // Optimizer and Objective Function
    private final MultivariateFunction tangentialHeadingTimeFunction;
    private final MultiStartMultivariateOptimizer optimizationPerformer;

    public TangentialHeadingSolver(
            double vMax, double mass, double muK, double c1, double c2, Point2D p0, Point2D p3,
            double thetaFinal, double boundaryTolerance, double submersibleDistanceTolerance, boolean isBlueAlliance,
            double thetaInitial, double angularVelocity, double robotWidth, double robotHeight) {
        this.v_max = vMax;
        this.mass = mass;
        this.mu_k = muK;
        this.c_1 = c1;
        this.c_2 = c2;
        this.p_0 = p0;
        this.p_3 = p3;
        this.theta_final = thetaFinal;
        this.isBlueAlliance = isBlueAlliance;
        this.theta_initial = thetaInitial;
        this.angularVelocity = angularVelocity;
        this.boundaryTolerance = boundaryTolerance;
        this.submersibleTolerance = submersibleDistanceTolerance;
        this.robotWidth = robotWidth;
        this.robotHeight = robotHeight;

        // Integrator
        this.integrator = new GillIntegrator(0.005);
        this.outputModel = new ContinuousOutputModel();

        // Define the optimization function
        this.tangentialHeadingTimeFunction = doubles -> {
            CubicBezierCurve bezierCurve = new CubicBezierCurve(
                    p_0,
                    new Point2D.Double(doubles[1], doubles[2]),
                    new Point2D.Double(doubles[3], doubles[4]),
                    p_3
            );
            return findT1(solveDifferentialEquation(doubles[0], bezierCurve), bezierCurve.getArcLength())
                    + findT2(doubles[0])
                    + bezierCurve.intersectionWeight(isBlueAlliance, doubles[0], boundaryTolerance, submersibleTolerance, robotWidth, robotHeight);
        };

        // Optimizer
        this.optimizationPerformer = new MultiStartMultivariateOptimizer(new BOBYQAOptimizer(10), 30, () -> {
            double[] output = new double[5];
            output[0] = Math.random() * (theta_final - theta_initial) + theta_initial;
            output[1] = Math.random() * (72 - 2 * boundaryTolerance) + boundaryTolerance;
            output[2] = Math.random() * (144 - 2 * boundaryTolerance) + boundaryTolerance;
            output[3] = Math.random() * (72 - 2 * boundaryTolerance) + boundaryTolerance;
            output[4] = Math.random() * (144 - 2 * boundaryTolerance) + boundaryTolerance;
            return output;
        });
    }

    private double[] performOptimization() {
        OptimizationData initialGuess = new InitialGuess(new double[]{
                theta_final, p_0.getX(), (p_0.getY() + p_3.getY()) / 2, p_0.getX(), p_3.getY()
        });
        OptimizationData bounds = new SimpleBounds(
                new double[]{0, boundaryTolerance, boundaryTolerance, boundaryTolerance, boundaryTolerance},
                new double[]{FastMath.PI * 2, 72 - boundaryTolerance, 144 - boundaryTolerance, 72 - boundaryTolerance, 144 - boundaryTolerance}
        );
        PointValuePair result = optimizationPerformer.optimize(
                initialGuess,
                new ObjectiveFunction(tangentialHeadingTimeFunction),
                GoalType.MINIMIZE,
                bounds
        );
        return result.getPoint();
    }

    private Point2D.Double[] solveDifferentialEquation(double theta, CubicBezierCurve bezierCurve) {
        TangentialHeadingDifferentialEquation diffeq = new TangentialHeadingDifferentialEquation(theta, v_max, mass, mu_k, c_1, c_2, bezierCurve);
        ExpandableStatefulODE expandableODE = new ExpandableStatefulODE(diffeq);
        expandableODE.setTime(0);
        expandableODE.setPrimaryState(new double[]{p_0.getX(), p_0.getY()});
        integrator.setMaxEvaluations(1000);
        integrator.addStepHandler(outputModel);
        integrator.integrate(expandableODE, 30.0);

        int n = 1000;
        Point2D.Double[] solutionPoints = new Point2D.Double[n];
        for (int i = 0; i < n; i++) {
            double t = (double) (i * 30) / 999.0;
            outputModel.setInterpolatedTime(t);
            double[] sol = outputModel.getInterpolatedState();
            solutionPoints[i] = new Point2D.Double(sol[0], sol[1]);
        }
        return solutionPoints;
    }

    private double findT1(Point2D.Double[] velocity, double L) {
        int minDistIndex = 0;
        double dist = Math.abs(velocity[minDistIndex].distance(0, 0) - L);

        for (int i = 1; i < velocity.length; i++) {
            double dist_i = Math.abs(velocity[i].distance(0, 0) - L);
            if (dist_i < dist) {
                minDistIndex = i;
                dist = dist_i;
            }
        }
        return (double) (minDistIndex * 30) / 999;
    }

    private double findT2(double theta) {
        return (FastMath.abs(theta_final - theta) + FastMath.abs(theta_initial - theta)) / angularVelocity;
    }

    public SolutionPoints getSolution() {
        double[] solution = performOptimization();
        double theta = solution[0];
        Point2D p_1 = new Point2D.Double(solution[1], solution[2]);
        Point2D p_2 = new Point2D.Double(solution[3], solution[4]);
        CubicBezierCurve path = new CubicBezierCurve(p_0, p_1, p_2, p_3);

        return new SolutionPoints(theta, path, findT1(solveDifferentialEquation(theta, path), path.getArcLength()), findT2(theta));
    }

    public static void main(String[] args) {
        TangentialHeadingSolver solver = new TangentialHeadingSolver(
                87, 16.09, 0.1, 10, 15,
                new Point2D.Double(10, 5),
                new Point2D.Double(10, 110),
                FastMath.PI / 2, 3, 3,
                true, 0, 6, 12.83, 15.75
        );
        SolutionPoints solution = solver.getSolution();
        System.out.println(solution);
    }
}