package MathUtil;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.solvers.LaguerreSolver;
import org.apache.commons.math3.complex.Complex;

import java.awt.geom.Point2D;

public class CubicBezierCurve {
    private Point2D P_0;
    private Point2D P_1;
    private Point2D P_2;
    private Point2D P_3;
    PolynomialFunction pathX;
    PolynomialFunction pathY;

    private LaguerreSolver solver = new LaguerreSolver();
    private SimpsonIntegrator integrator = new SimpsonIntegrator();

    private double arcLength = -1;

    public CubicBezierCurve(Point2D P_0, Point2D P_1, Point2D P_2, Point2D P_3) {
        this.P_0 = P_0;
        this.P_1 = P_1;
        this.P_2 = P_2;
        this.P_3 = P_3;

        double cubeTermX = P_3.getX() - 3*P_2.getX() + 3*P_1.getX() - P_0.getX();
        double cubeTermY = P_3.getY() - 3*P_2.getY() + 3*P_1.getY() - P_0.getY();
        double squareTermX = 3*P_2.getX() - 6*P_1.getX() + 3*P_0.getX();
        double squareTermY = 3*P_2.getY() - 6*P_1.getY() + 3*P_0.getY();
        double proportionalTermX = 3*P_1.getX()-3*P_0.getX();
        double proportionalTermY = 3*P_1.getY()-3*P_0.getY();
        double constantTermX = P_0.getX();
        double constantTermY = P_0.getY();

        pathX = new PolynomialFunction(new double[] {cubeTermX, squareTermX, proportionalTermX, constantTermX});
        pathY = new PolynomialFunction(new double[] {cubeTermY, squareTermY, proportionalTermY, constantTermY});
    }

    public Point2D.Double compute(double t) {
        return new Point2D.Double(pathX.value(t), pathY.value(t));
    }

    public Point2D.Double computePathPrimes(double t) {
        return new Point2D.Double(pathX.polynomialDerivative().value(t), pathY.polynomialDerivative().value(t));
    }

    public double pathInversion(Point2D point) {
        Complex[] t_x = solver.solveAllComplex(pathX.getCoefficients(), 0.5, 300);
        Complex[] t_y = solver.solveAllComplex(pathY.getCoefficients(), 0.5, 300);

        double root = 0.5;

        for (Complex root_x : t_x) {
            for (Complex root_y : t_y) {
                if (root_x == root_y && root_x.getImaginary() == 0 && root_x.getReal() <= 1 && root_x.getReal() >= 0) {
                    root = root_x.getReal();
                }
            }
        }

        return root;
    }

    public double getArcLength() {
        if (arcLength == -1) {
            UnivariateFunction arcFunction = t -> computePathPrimes(t).distance(0,0);
            arcLength = integrator.integrate(200, arcFunction, 0, 1);
        }

        return arcLength;
    }

    public Point2D getP_0() {
        return P_0;
    }

    public Point2D getP_1() {
        return P_1;
    }

    public Point2D getP_2() {
        return P_2;
    }

    public Point2D getP_3() {
        return P_3;
    }
}