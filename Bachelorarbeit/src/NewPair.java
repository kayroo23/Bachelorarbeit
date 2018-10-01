class NewPair extends Point{
    private double[] attribute = new double[2];

    NewPair(int x, int y){
        super(2);
        attribute[0] = x;
        attribute[1] = y;
    }

    double[] getAttribute(){
        return attribute;
    }

    double getFirst(){
        return attribute[0];
    }

    double getSecond(){
        return attribute[1];
    }
}
