package progistar.pXg.data;

public class Cigar {
	public int		markerSize;
	public char		operation;
	public String	nucleotides = ""; // default

	// relative position to the start site of NGS-read
	public int[]	relativePositions;

	// see Regional Characters in Constants Class
	// char[relativePositionSize][txdIndiciesSize]
	// each relative position has regional information corresponding to the transcript.
	public char[][]	annotations;

	public Cigar (int markerSize, char operation) {
		this.markerSize = markerSize;
		this.operation = operation;
	}
}
