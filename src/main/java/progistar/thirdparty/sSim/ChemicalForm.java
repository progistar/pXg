package progistar.thirdparty.sSim;


public class ChemicalForm {

	private int S = 0;
	private int N = 0;
	private int O = 0;
	private int H = 0;
	private int C = 0;
	private int P = 0;
	private int Na = 0;
	private int Se = 0;
	private int Fe = 0;
	private int K = 0;
	private int I = 0;

	public ChemicalForm () {

	}

	public ChemicalForm (int C, int H, int O, int N, int S) {
		this.H = H;
		this.C = C;
		this.N = N;
		this.O = O;
		this.S = S;
	}

	public ChemicalForm (int C, int H, int O, int N, int S, int P) {
		this.H = H;
		this.C = C;
		this.N = N;
		this.O = O;
		this.S = S;
		this.P = P;
	}

	public ChemicalForm (int C, int H, int O, int N, int S, int P, int Se) {
		this.H = H;
		this.C = C;
		this.N = N;
		this.O = O;
		this.S = S;
		this.P = P;
		this.Se = Se;
	}

	public void setS (int count) {
		this.S = count;
	}

	public void setN (int count) {
		this.N = count;
	}

	public void setO (int count) {
		this.O = count;
	}

	public void setH (int count) {
		this.H = count;
	}

	public void setC (int count) {
		this.C = count;
	}

	public void setP (int count) {
		this.P = count;
	}

	public void setNa (int count) {
		this.Na = count;
	}

	public void setSe (int count) {
		this.Se = count;
	}

	public void setFe (int count) {
		this.Fe = count;
	}

	public void setI (int count) {
		this.I = count;
	}

	public void setK (int count) {
		this.K = count;
	}

	public int getC () {
		return this.C;
	}

	public int getH () {
		return this.H;
	}

	public int getO () {
		return this.O;
	}

	public int getN () {
		return this.N;
	}

	public int getS () {
		return this.S;
	}

	public int getP () {
		return this.P;
	}

	public int getNa () {
		return this.Na;
	}

	public int getSe () {
		return this.Se;
	}

	public int getFe() {
		return this.Fe;
	}

	public int getK () {
		return this.K;
	}

	public int getI () {
		return this.I;
	}

	public double getMass () {
		double mass = ProteomeConstants.Carbon * this.C +
				ProteomeConstants.Hydrogen * this.H +
				ProteomeConstants.Oxygen * this.O +
				ProteomeConstants.Nitrogen * this.N +
				ProteomeConstants.Sulfur * this.S +
				ProteomeConstants.Phosphorus * this.P +
				ProteomeConstants.Sodium * this.Na +
				ProteomeConstants.Selenium * this.Se +
				ProteomeConstants.Ferrum * this.Fe +
				ProteomeConstants.Potassium * this.K +
				ProteomeConstants.Iodine * this.I;
		return mass;
	}

	@Override
	public String toString() {
		StringBuilder formula = new StringBuilder();
		if(this.C > 0) {
			formula.append("C").append("(").append(this.C).append(")");
		}
		if(this.H > 0) {
			formula.append("H").append("(").append(this.H).append(")");
		}
		if(this.O > 0) {
			formula.append("O").append("(").append(this.O).append(")");
		}
		if(this.N > 0) {
			formula.append("N").append("(").append(this.N).append(")");
		}
		if(this.S > 0) {
			formula.append("S").append("(").append(this.S).append(")");
		}
		if(this.P > 0) {
			formula.append("P").append("(").append(this.P).append(")");
		}
		if(this.Na > 0) {
			formula.append("Na").append("(").append(this.Na).append(")");
		}
		if(this.Se > 0) {
			formula.append("Se").append("(").append(this.Se).append(")");
		}
		if(this.Fe > 0) {
			formula.append("Fe").append("(").append(this.Fe).append(")");
		}
		if(this.K > 0) {
			formula.append("K").append("(").append(this.K).append(")");
		}
		if(this.I > 0) {
			formula.append("I").append("(").append(this.I).append(")");
		}
		return formula.toString();
	}


}