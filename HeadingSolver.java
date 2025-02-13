import ConstantHeadingInterpolation.ConstantHeadingDifferentialEquation;
import LinearHeadingInterpolation.LinearHeadingDifferentialEquation;
import MathUtil.CubicBezierCurve;
import MathUtil.SolutionPoints;
import TangentialHeadingInterpolation.TangentialHeadingDifferentialEquation;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.ExpandableStatefulODE;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.nonstiff.GillIntegrator;
import org.apache.commons.math3.util.FastMath;

import java.awt.geom.Point2D;

public class HeadingSolver {
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
    private final double boundaryTolerance;
    private final double submersibleTolerance;
    private final double robotWidth;
    private final double robotHeight;
    private final double f_heading;
    private final double d_heading;
    private final double p_heading;

    private GillIntegrator integrator = new GillIntegrator(0.005);
    private ContinuousOutputModel outputModel = new ContinuousOutputModel();

    public enum HeadingState {
        CONSTANT,
        LINEAR,
        TANGENTIAL
    }
    private HeadingState headingState = HeadingState.CONSTANT;

    private MultivariateFunction timeFunction;

    public HeadingSolver(double vMax, double mass, double muK, double c1, double c2, Point2D p0, Point2D p3, double thetaFinal, boolean isBlueAlliance, double thetaInitial, double angularVelocity, double boundaryTolerance, double submersibleTolerance, double robotWidth, double robotHeight, double fHeading, double dHeading, double pHeading) {
        v_max = vMax;
        this.mass = mass;
        mu_k = muK;
        c_1 = c1;
        c_2 = c2;
        p_0 = p0;
        p_3 = p3;
        theta_final = thetaFinal;
        this.isBlueAlliance = isBlueAlliance;
        theta_initial = thetaInitial;
        this.angularVelocity = angularVelocity;
        this.boundaryTolerance = boundaryTolerance;
        this.submersibleTolerance = submersibleTolerance;
        this.robotWidth = robotWidth;
        this.robotHeight = robotHeight;
        f_heading = fHeading;
        d_heading = dHeading;
        p_heading = pHeading;
        timeFunction = doubles -> {
            CubicBezierCurve bezierCurve = new CubicBezierCurve(p_0, new Point2D.Double(doubles[1], doubles[2]), new Point2D.Double(doubles[3], doubles[4]),p_3);
            return findT1(solveDifferentialEquation(doubles[0], bezierCurve),bezierCurve.getArcLength()) + findT2(doubles[0]);
        };
    }

    public void setHeadingState(HeadingState headingState) {
        this.headingState = headingState;
    }

    public HeadingState getHeadingState() {
        return headingState;
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

    private double findT2(double deltaHeading) {
        return FastMath.abs(theta_final - (deltaHeading + theta_initial))/angularVelocity;
    }

    public SolutionPoints getSolution(HeadingState headingState) {
        this.headingState = headingState;
        double[] solution = performOptimization();
        double theta = solution[0];
        Point2D p_1 = new Point2D.Double(solution[1], solution[2]);
        Point2D p_2 = new Point2D.Double(solution[3], solution[4]);
        CubicBezierCurve path = new CubicBezierCurve(p_0, p_1, p_2, p_3);

        return new SolutionPoints(theta, path, findT1(solveDifferentialEquation(theta, path), path.getArcLength()), findT2(theta));
    }

    private double[] performOptimization() {
        //TODO: Complete this method
        return null;
    }

    private Point2D.Double[] solveDifferentialEquation(double thetaClassifier, CubicBezierCurve bezierCurve) {
        ExpandableStatefulODE expandableODE = null;
        if (headingState == HeadingState.CONSTANT) {
            ConstantHeadingDifferentialEquation diffeq = new ConstantHeadingDifferentialEquation(thetaClassifier, v_max, mass, mu_k, c_1, bezierCurve);
            expandableODE = new ExpandableStatefulODE(diffeq);
        } else if (headingState == HeadingState.LINEAR) {
            LinearHeadingDifferentialEquation diffeq = new LinearHeadingDifferentialEquation(thetaClassifier, theta_initial, v_max, mass, mu_k, c_1, c_2, f_heading, d_heading, p_heading, angularVelocity, bezierCurve);
            expandableODE = new ExpandableStatefulODE(diffeq);
        } else if (headingState == HeadingState.TANGENTIAL) {
            TangentialHeadingDifferentialEquation diffeq = new TangentialHeadingDifferentialEquation(v_max, mass, mu_k, c_2, bezierCurve, angularVelocity, f_heading, d_heading, p_heading);
            expandableODE = new ExpandableStatefulODE(diffeq);
        }
        
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
}
