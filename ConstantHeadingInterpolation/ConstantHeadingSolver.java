package ConstantHeadingInterpolation;

import MathUtil.CubicBezierCurve;
import MathUtil.SolutionPoints;
import MathUtil.SquareCoordinates;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.ExpandableStatefulODE;
import org.apache.commons.math3.ode.nonstiff.GillIntegrator;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.OptimizationData;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateFunctionMappingAdapter;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import org.apache.commons.math3.util.FastMath;

import java.awt.geom.Point2D;

public class ConstantHeadingSolver {
    private final double v_max; //inches per second
    private final double mass; //kg
    private final double mu_k;
    private final double c_1;
    private final double c_2;
    private final Point2D p_0;
    private final Point2D p_3;
    private final double theta_final;
    private final boolean isBlueAlliance;
    private final double theta_initial;
    private final double angularVelocity;

    private GillIntegrator integrator = new GillIntegrator(0.005);
    private ContinuousOutputModel outputModel = new ContinuousOutputModel();

    private MultivariateFunction constantHeadingTimeFunction = new MultivariateFunction() {
        @Override
        public double value(double[] doubles) {
            CubicBezierCurve bezierCurve = new CubicBezierCurve(p_0, new SquareCoordinates(doubles[1], doubles[2]).getCartesian(), new SquareCoordinates(doubles[3],doubles[4]).getCartesian(),p_3);
           return findT1(solveDifferentialEquation(doubles[0], bezierCurve),bezierCurve.getArcLength()) + findT2(doubles[0]);

        }
    };

    //doubles[0] = theta
    //doubles[1, 2] = square coordinates p_1
    //doubles[3, 4] = square coordinates p_2

    private MultivariateFunction time;

    private BOBYQAOptimizer optimizer = new BOBYQAOptimizer(10000);

    public ConstantHeadingSolver(double vMax, double mass, double muK, double c1, double c2, Point2D p0, Point2D p3, double thetaT2, double boundaryTolerance, double submersibleDistanceTolerance, boolean isBlueAlliance, double thetaInitial, double angularVelocity) {
        v_max = vMax;
        this.mass = mass;
        mu_k = muK;
        c_1 = c1;
        c_2 = c2;
        p_0 = p0;
        p_3 = p3;
        theta_final = thetaT2;
        this.isBlueAlliance = isBlueAlliance;
        theta_initial = thetaInitial;
        this.angularVelocity = angularVelocity;

        if (isBlueAlliance) {
            time = new MultivariateFunctionMappingAdapter(constantHeadingTimeFunction, new double[]{0, submersibleDistanceTolerance, 0, submersibleDistanceTolerance, 0}, new double[]{2 * FastMath.PI, 72 - boundaryTolerance, FastMath.PI, 72 - boundaryTolerance, FastMath.PI});
        } else {
            time = new MultivariateFunctionMappingAdapter(constantHeadingTimeFunction, new double[]{0, submersibleDistanceTolerance, Math.PI, submersibleDistanceTolerance, Math.PI}, new double[]{2 * FastMath.PI, 72 - boundaryTolerance, 2 * FastMath.PI, 72 - boundaryTolerance, 2 * FastMath.PI});
        }
    }

    private double[] performOptimization() {
        OptimizationData optimizationData = new InitialGuess(new double[] {theta_final, p_0.getX(), (p_0.getY() + p_3.getY())/2, p_0.getX(), p_3.getY()});
        PointValuePair result = optimizer.optimize(optimizationData, new ObjectiveFunction(time), GoalType.MINIMIZE);
        return result.getPoint();
    }

    private Point2D.Double[] solveDifferentialEquation(double theta, CubicBezierCurve bezierCurve) {
        ConstantHeadingDifferentialEquation diffeq = new ConstantHeadingDifferentialEquation(theta, v_max, mass, mu_k, c_1, c_2, bezierCurve);
        ExpandableStatefulODE expandableODE = new ExpandableStatefulODE(diffeq);
        expandableODE.setTime(0);
        double[] primaryState = {p_0.getX(),p_0.getY()};
        expandableODE.setPrimaryState(primaryState);
        integrator.setMaxEvaluations(1000);
        integrator.integrate(expandableODE, 30.0); //if you change lines 62 or 63 then also change the denominator or numerator on line 69 respectively

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
        int min_dist_index = 0;
        double dist = Math.abs(velocity[min_dist_index].distance(0,0) - L);

        for (int i = 0; i < velocity.length; i++) {
            Point2D point = velocity[i];
            double dist_i = Math.abs(point.distance(0,0) - L);
            if (dist_i < dist) {
                min_dist_index = i;
                dist = dist_i;
            }
        }

        return (double) (min_dist_index * 30)/999;
    }

    private double findT2(double theta) {
        return (FastMath.abs(theta_final - theta) + FastMath.abs(theta_initial - theta))/angularVelocity;
    }

   public SolutionPoints getSolution() {
        double[] solution = performOptimization();
        double theta = solution[0];
        SquareCoordinates p_1 = new SquareCoordinates(solution[1], solution[2]);
        SquareCoordinates p_2 = new SquareCoordinates(solution[3], solution[4]);
        CubicBezierCurve path = new CubicBezierCurve(p_0, p_1.getCartesian(), p_2.getCartesian(), p_3);

        return new SolutionPoints(theta, path, findT1(solveDifferentialEquation(theta, path), path.getArcLength()), findT2(theta));
   }
}
