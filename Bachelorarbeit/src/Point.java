import java.util.ArrayList;

public class Point {

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



}
