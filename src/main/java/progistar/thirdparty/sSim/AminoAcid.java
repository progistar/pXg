package progistar.thirdparty.sSim;

import java.util.ArrayList;


public class AminoAcid {

	private static final AminoAcid [] AMINO_ACID_TABLE =
		{
			new AminoAcid('A',	71.037114,	"alanine", new ChemicalForm(3, 5, 1, 1, 0)),
			null,
			new AminoAcid('C',	103.00919,	"cysteine", new ChemicalForm(3, 5, 1, 1, 1)),
			new AminoAcid('D',	115.02694,	"aspartate", new ChemicalForm(4, 5, 3, 1, 0)),
			new AminoAcid('E',	129.04259,	"glutamate", new ChemicalForm(5, 7, 3, 1, 0)),
			new AminoAcid('F',	147.06841,	"phenylalanine", new ChemicalForm(9, 9, 1, 1, 0)),
			new AminoAcid('G',	57.021464,	"glycine", new ChemicalForm(2, 3, 1, 1, 0)),
			new AminoAcid('H',	137.05891,	"histidine", new ChemicalForm(6, 7, 1, 3, 0)),
			new AminoAcid('I',	113.08406,	"isoleucine", new ChemicalForm(6, 11, 1, 1, 0)),
			null,
			new AminoAcid('K',	128.09496,	"lysine", new ChemicalForm(6, 12, 1, 2, 0)),
			new AminoAcid('L',	113.08406,	"leucine", new ChemicalForm(6, 11, 1, 1, 0)),
			new AminoAcid('M',	131.04048,	"methionine", new ChemicalForm(5, 9, 1, 1, 1)),
			new AminoAcid('N',	114.04293,	"asparagine", new ChemicalForm(4, 6, 2, 2, 0)),
			null,
			new AminoAcid('P',	97.052764,	"proline", new ChemicalForm(5, 7, 1, 1, 0)),
			new AminoAcid('Q',	128.05858,	"glutamine", new ChemicalForm(5, 8, 2, 2, 0)),
			new AminoAcid('R',	156.10111,	"arginine", new ChemicalForm(6, 12, 1, 4, 0)),
			new AminoAcid('S',	87.032029,	"serine", new ChemicalForm(3, 5, 2, 1, 0)),
			new AminoAcid('T',	101.04768,	"threonine", new ChemicalForm(4, 7, 2, 1, 0)),
			new AminoAcid('U',	151.00919,	"selenocysteine", new ChemicalForm(3, 5, 1, 1, 0, 0, 1)),
			new AminoAcid('V',	99.068414,	"valine", new ChemicalForm(5, 9, 1, 1, 0)),
			new AminoAcid('W',	186.07931,	"tryptophan", new ChemicalForm(11, 10, 1, 2, 0)),
			null,
			new AminoAcid('Y',	163.06333,	"tyrosine", new ChemicalForm(9, 9, 2, 1, 0)),
			null,
		};


	private char aminoAcid = 'X';
	private double mass = .0;
	private String fullName = null;
	private ChemicalForm chemicalForm = null;

	public AminoAcid (char aminoAcid, double mass, String fullName, ChemicalForm chemicalForm) {
		this.aminoAcid = aminoAcid;
		this.mass = mass;
		this.fullName = fullName;
		this.chemicalForm = chemicalForm;
	}

	public ChemicalForm getChemicalForm () {
		return this.chemicalForm;
	}

	public double getResidualMass () {
		return this.mass;
	}

	public char getAminoAcid () {
		return this.aminoAcid;
	}

	public String getFullName () {
		return this.fullName;
	}

	public static AminoAcid getAminoAcid (char aminoAcid) {
		if(Character.isLowerCase(aminoAcid)) {
			aminoAcid = Character.toUpperCase(aminoAcid);
		}
		return AMINO_ACID_TABLE[aminoAcid - 'A'];
	}

	/**
	 * tolerance must be dalton. <br>
	 * If there is no matched amino acid, then it return zero size arraylist. <br>
	 * @param mass
	 * @param charge
	 * @param tolerance
	 * @return
	 */
	public static ArrayList<AminoAcid> getAminoAcid (double mass, double tolerance) {
		ArrayList<AminoAcid> aminoAcids = new ArrayList<>();

		for (AminoAcid element : AminoAcid.AMINO_ACID_TABLE) {
			if(element != null) {
				if(Math.abs(mass - element.getResidualMass()) < tolerance) {
					aminoAcids.add(element);
				}
			}
		}

		return aminoAcids;
	}

	/**
	 * tolerance must be dalton. <br>
	 * If there is no matched amino acid, then it returns null. <br>
	 *
	 * @param mass
	 * @param charge
	 * @param tolerance
	 * @param aminoAcid
	 * @return
	 */
	public static AminoAcid getSpecificAminoAcid (double mass, double tolerance, char aminoAcid) {
		AminoAcid aminoAcid_ = null;
		int aminoIndex = aminoAcid - 'A';
		if(AminoAcid.AMINO_ACID_TABLE[aminoIndex] != null) {
			if(Math.abs(mass - AminoAcid.AMINO_ACID_TABLE[aminoIndex].getResidualMass()) < tolerance) {
				aminoAcid_ = AminoAcid.AMINO_ACID_TABLE[aminoIndex];
			}
		}

		return aminoAcid_;
	}
}