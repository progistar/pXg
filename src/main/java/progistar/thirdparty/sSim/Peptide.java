package progistar.thirdparty.sSim;

import java.util.ArrayList;
import java.util.Arrays;


/**
 * Peptide class covers almost possible characteristics, which can be occurred in peptide. <br>
 * For example, several PTMs and variants (SNV, INDEL) can be considered. <br>
 * 
 * 
 * @author progi
 *
 */
public class Peptide {
	private String sequence = null;
	private String preAA = null;
	private String postAA = null;
	private ArrayList<Modification>[] modifications = null;
	
	public Peptide (String sequence, String preAA, String postAA) {
		this.sequence = sequence;
		this.preAA = preAA;
		this.postAA = postAA;
		this.modifications = new ArrayList[sequence.length()+1];
	}

	public Peptide getDecoyPeptide () {
		Peptide decoyPeptide = new Peptide(new StringBuilder(this.sequence).reverse().toString(), this.postAA, this.preAA);
		
		//reverse modifications (nterm is fixed)
		for(int i=1; i<sequence.length()+1; i++) decoyPeptide.modifications[i] = this.modifications[sequence.length()-i+1];
		
		return decoyPeptide;
	}
	
	/**
	 * ion: The available ion list is represented in ProteomeConstants <br>
	 * 
	 * @param ion
	 * @param charge
	 * @param withMods
	 * @return
	 */
	public double[] getTheoreticalLadder (double ion, double charge, boolean withMods) {
		int sequenceLength = sequence.length();
		double[] ladder = new double[sequenceLength];
		
		double nTermModi = 0;
		ArrayList<Modification> nTermMods = this.getModification(0);
		if(withMods && nTermMods != null) for(Modification mod : nTermMods) nTermModi += mod.getModMass();
		
		for(int i=0; i<sequenceLength; i++) {
			ladder[i] = AminoAcid.getAminoAcid(sequence.charAt(i)).getResidualMass();
			ArrayList<Modification> mods = this.getModification(i+1);
			if(withMods && mods != null) for(Modification mod : mods) {
				ladder[i] += mod.getModMass();
			}
		}
		
		// x, y, z ions
		if(ion == ProteomeConstants.X_ION || ion == ProteomeConstants.Y_ION || ion == ProteomeConstants.Z_ION) {
			// accumulated ions
			for(int i=sequenceLength-2; i>=0; i--) ladder[i] += ladder[i+1];
		}
		// a, b, c ions
		else if (ion == ProteomeConstants.A_ION || ion == ProteomeConstants.B_ION || ion == ProteomeConstants.C_ION){
			// accumulated ions
			ladder[0] += nTermModi;
			for(int i=1; i<sequenceLength; i++) ladder[i] += ladder[i-1];
		} else {
			// unsupported ions
			System.err.println("Unsupported IONs: "+ion);
			return null;
		}
		
		// charge and ion
		for(int i=0; i<sequenceLength; i++) ladder[i] = (ladder[i] + (charge-1)*ProteomeConstants.Proton + ion) / charge;
		
		Arrays.sort(ladder);
		
		double[] refinedLadder = new double[ladder.length-1];
		for(int i=0; i<refinedLadder.length; i++) {
			refinedLadder[i] = ladder[i];
		}
		
		return refinedLadder;
	}
	
	/**
	 * Position info <br>
	 * N-term: 0 <br>
	 * Others: one-based <br>
	 * 
	 * @param mod
	 * @param position
	 */
	public void addModification (Modification mod, int position) {
		if(this.modifications[position] == null) this.modifications[position] = new ArrayList<Modification>();
		this.modifications[position].add(mod);
	}
	
	/**
	 * Position info <br>
	 * N-term: 0 <br>
	 * Others: one-based<br>
	 * @param position
	 * @return
	 */
	public ArrayList<Modification> getModification (int position) {
		return this.modifications[position];
	}
	
	public String getSequence (boolean withMods) {
		if(withMods) {
			StringBuilder modifiedPeptide = new StringBuilder(this.sequence).append("[");
			boolean isModified = false;
			for(int pos = 0; pos < this.modifications.length; pos++) {
				ArrayList<Modification> mods = this.getModification(pos);
				if(mods != null) for(Modification mod : mods) {
					if(isModified) modifiedPeptide.append(","); 
					modifiedPeptide.append(pos+":"+mod.modName);
					isModified = true;
				}
			}
			
			return modifiedPeptide.append("]").toString();
			
		} else return this.sequence;
	}
	
	public void setSequence (String sequence) {
		this.sequence = sequence;
	}
	
	public String getPreAA () {
		return preAA;
	}
	
	public void setPreAA (String preAA) {
		this.preAA = preAA;
	}
	
	public String getPostAA () {
		return this.postAA;
	}
	
	public void setPostAA (String postAA) {
		this.postAA = postAA;
	}

	/**
	 * mass of peptide with H20.<br>
	 * 
	 * @param withMods
	 * @return
	 */
	public double getTheoreticalMass (boolean withMods) {
		double theoreticalMass = ProteomeConstants.H2O;
		for(int i=0; i<sequence.length(); i++) theoreticalMass += AminoAcid.getAminoAcid(sequence.charAt(i)).getResidualMass();
		if(withMods) {
			for(ArrayList<Modification> mods : this.modifications) {
				if(mods != null) for(Modification mod : mods) theoreticalMass += Double.parseDouble(mod.modMass);
			}
		}
		
		return theoreticalMass;
	}
	
	/**
	 * chemical formula of peptide with H2O.<br>
	 * 
	 * @param withMods
	 * @return
	 */
	public String getChemicalFormString (boolean withMods) {
		ChemicalForm thisChemicalForm = new ChemicalForm();
		thisChemicalForm.setH(2); thisChemicalForm.setO(1); // set H2O
		int length = this.sequence.length();
		for(int i=0; i<length; i++) {
			ChemicalForm AAChem = AminoAcid.getAminoAcid(sequence.charAt(i)).getChemicalForm();
			thisChemicalForm.setC(thisChemicalForm.getC() + AAChem.getC());
			thisChemicalForm.setH(thisChemicalForm.getH() + AAChem.getH());
			thisChemicalForm.setO(thisChemicalForm.getO() + AAChem.getO());
			thisChemicalForm.setN(thisChemicalForm.getN() + AAChem.getN());
			thisChemicalForm.setS(thisChemicalForm.getS() + AAChem.getS());
		}
		
		if(withMods) {
			for(ArrayList<Modification> mods : this.modifications) {
				if(mods != null) {
					for(Modification mod : mods) {
						ChemicalForm modChem = mod.getChemicalForm();
						if(modChem == null) {
							System.err.println("Fail to get chemical string. There is unknown modification.");
							return null;
						} else {
							thisChemicalForm.setC(thisChemicalForm.getC() + modChem.getC());
							thisChemicalForm.setH(thisChemicalForm.getH() + modChem.getH());
							thisChemicalForm.setO(thisChemicalForm.getO() + modChem.getO());
							thisChemicalForm.setN(thisChemicalForm.getN() + modChem.getN());
							thisChemicalForm.setS(thisChemicalForm.getS() + modChem.getS());
							thisChemicalForm.setP(thisChemicalForm.getP() + modChem.getP());
							thisChemicalForm.setNa(thisChemicalForm.getNa() + modChem.getNa());
							thisChemicalForm.setSe(thisChemicalForm.getSe() + modChem.getSe());
							thisChemicalForm.setFe(thisChemicalForm.getFe() + modChem.getFe());
							thisChemicalForm.setK(thisChemicalForm.getK() + modChem.getK());
							
						}
					}
				}
			}
		}
		
		return thisChemicalForm.toString();
	}
	
}