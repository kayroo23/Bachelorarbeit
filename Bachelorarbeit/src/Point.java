import org.apache.commons.math3.ml.clustering.Clusterable;

public class Point implements Clusterable {

    private int numberOfAttributes;
    private double[] attributes;

    public Point(int numberOfAttributes){
        this.numberOfAttributes = numberOfAttributes;
        attributes = new double[this.numberOfAttributes];
    }

    public Point(int numberOfAttributes, double[] att){
        this.numberOfAttributes = numberOfAttributes;
        attributes = new double[this.numberOfAttributes];
        this.attributes = att;
    }

    int getNumberOfAttributes(){
        return numberOfAttributes;
    }

    double[] getAttributes(){
        return attributes;
    }

    public double[] getPoint(){
        return attributes;
    }

}
