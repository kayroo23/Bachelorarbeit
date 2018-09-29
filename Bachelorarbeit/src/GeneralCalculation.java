import org.apache.commons.math3.linear.RealMatrix;

import java.util.List;

public abstract class GeneralCalculation {

    protected abstract double[][] calculateMatrix(List<List<Point>> newClusters, int numberOfPointsPerCluster);
    protected abstract List<Integer> calculateBestAttributes(int numberOfShownAttributes, RealMatrix spearmanMatrix);
    protected abstract void calculateResults(List<List<Point>> newClusters, List<Integer> bestAttributes);

}
