package progistar.thirdparty.sSim;


import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Iterator;

public class Modification {

	public static final String ANY_NTERM = "Any-Nterm";
	public static final String NTERM = "Nterm";

	protected String modName = null;
	protected String modMass = null;
	protected String modSite = null;
	protected ChemicalForm chemicalForm = null;
	protected boolean isFixed = false;

	private static Hashtable<String, ArrayList<Modification>> preDefinedModificationList = null;
	private static Hashtable<String, String> massByModName = null;

	public Modification (String modName, String modSite, String modMass, boolean isFixed) {
		this.modName = modName;
		this.modSite = modSite;
		this.modMass = modMass;
		this.isFixed = isFixed;
	}

	public Modification (String modName, String modSite, String modMass, boolean isFixed, ChemicalForm chemicalForm) {
		this.modName = modName;
		this.modSite = modSite;
		this.modMass = modMass;
		this.isFixed = isFixed;
		this.chemicalForm = chemicalForm;
	}

	/**
	 * Input example1: "K+229.163"<br>
	 * Input example2: "Modification.ANY_NTERM+229.163" <br>
	 * Input example3: "Modification.NTERM+229.163K" <br>
	 *
	 *
	 * @param modIdentifier
	 * @return
	 */
	public static ArrayList<Modification> getModificationsByModIdentifier (String modIdentifier) {
		if(preDefinedModificationList == null) {
			preDefinedModificationList = new Hashtable<>();
			massByModName = new Hashtable<>();
			initModification();
		}

		return preDefinedModificationList.get(modIdentifier);
	}

	public static String getModificationMassByModName (String modName) {
		if(preDefinedModificationList == null) {
			preDefinedModificationList = new Hashtable<>();
			massByModName = new Hashtable<>();
			initModification();
		}

		return massByModName.get(modName);
	}

	private static void initModification () {
		ChemicalForm chemicalForm = null;

		Modification TMT_NTERM = new Modification("Nterm-TMT", ANY_NTERM, "+229.163", true, new ChemicalForm(12, 20, 2, 2, 0)); // the actual form is H(20) C(8) 13C(4) N 15N O(2)
		Modification TMT_K = new Modification("Lysine-TMT", "K", "+229.163", true, new ChemicalForm(12, 20, 2, 2, 0));
		Modification CARBAM_C = new Modification("Cysteine-Carbamidomethly", "C", "+57.021", true, new ChemicalForm(2, 3, 1, 1, 0));
		Modification CARBAM_NTERM = new Modification("Nterm-Carbamidomethly", ANY_NTERM, "+57.021", false, new ChemicalForm(2, 3, 1, 1, 0));
		Modification CARBAM_K = new Modification("Lysine-Carbamidomethly", "K", "+57.021", false, new ChemicalForm(2, 3, 1, 1, 0));
		Modification CARBAM_H = new Modification("Histidine-Carbamidomethly", "H", "+57.021", false, new ChemicalForm(2, 3, 1, 1, 0));

		Modification CARBAMDTT_C = new Modification("Cysteine-CarbamidomethlyDTT", "C", "+209.018", false, new ChemicalForm(6, 11, 3, 1, 2));

		Modification CARBAMYL_K = new Modification("Lysine-Carbamyl", "K", "+43.006", false, new ChemicalForm(1, 1, 1, 1, 0));
		Modification CARBAMYL_NTERM = new Modification("Nterm-Carbamyl", ANY_NTERM, "+43.006", false, new ChemicalForm(1, 1, 1, 1, 0));

		chemicalForm = new ChemicalForm(0, -1, 0, 0, 0);
		chemicalForm.setNa(1);
		Modification CATION_NA_E = new Modification("GlutamincAcid-Cation:Na", "E", "+21.982", false, chemicalForm);
		Modification CATION_NA_D = new Modification("AsparticAcid-Cation:Na", "D", "+21.982", false, chemicalForm);
		chemicalForm = new ChemicalForm(0, -2, 0, 0, 0);
		chemicalForm.setFe(1);
		Modification CATION_FE_II_E = new Modification("GlutamincAcid-Cation:Fe[II]", "E", "+53.919", false, chemicalForm);
		Modification CATION_FE_II_D = new Modification("AsparticAcid-Cation:Fe[II]", "D", "+53.919", false, chemicalForm);
		chemicalForm = new ChemicalForm(0, -1, 0, 0, 0);
		chemicalForm.setK(1);
		Modification CATION_K_E = new Modification("GlutamincAcid-Cation:K", "E", "+37.956", false, chemicalForm);
		Modification CATION_K_D = new Modification("AsparticAcid-Cation:K", "D", "+37.956", false, chemicalForm);

		Modification OXIDATION_M = new Modification("Methionine-Oxidation", "M", "+15.995", false, new ChemicalForm(0, 0, 1, 0, 0));
		Modification OXIDATION_P = new Modification("Proline-Oxidation", "P", "+15.995", false, new ChemicalForm(0, 0, 1, 0, 0));
		Modification OXIDATION_W = new Modification("Tryptophan-Oxidation", "W", "+15.995", false, new ChemicalForm(0, 0, 1, 0, 0));
		Modification OXIDATION_C = new Modification("Cysteine-Oxidation", "C", "+15.995", false, new ChemicalForm(0, 0, 1, 0, 0));
		Modification DIOXIDATION_M = new Modification("Methionine-Dioxidation", "M", "+31.990", false, new ChemicalForm(0, 0, 2, 0, 0));
		Modification DIOXIDATION_C = new Modification("Cysteine-Dioxidation", "C", "+31.990", false, new ChemicalForm(0, 0, 2, 0, 0));
		Modification DIOXIDATION_W = new Modification("Tryptophan-Dioxidation", "W", "+31.990", false, new ChemicalForm(0, 0, 2, 0, 0));
		Modification TRIOXIDATION_C = new Modification("Cysteine-Trioxidation", "C", "+47.985", false, new ChemicalForm(0, 0, 3, 0, 0));

		Modification METHYL_K = new Modification("Lysine-Methylation", "K", "+14.016", false, new ChemicalForm(1, 2, 0, 0, 0));
		Modification METHYL_R = new Modification("Arginine-Methylation", "R", "+14.016", false, new ChemicalForm(1, 2, 0, 0, 0));
		Modification METHYL_E = new Modification("GlutamicAcid-Methylation", "E", "+14.016", false, new ChemicalForm(1, 2, 0, 0, 0));
		Modification METHYL_D = new Modification("AsparticAcid-Methylation", "D", "+14.016", false, new ChemicalForm(1, 2, 0, 0, 0));
		Modification METHYL_H = new Modification("Histidine-Methylation", "H", "+14.016", false, new ChemicalForm(1, 2, 0, 0, 0));
		Modification DIMETHYL_K = new Modification("Lysine-Dimethylation", "K", "+28.031", false, new ChemicalForm(2, 4, 0, 0, 0));
		Modification DIMETHYL_R = new Modification("Arginine-Dimethylation", "R", "+28.031", false, new ChemicalForm(2, 4, 0, 0, 0));
		Modification TRIMETHYL_K = new Modification("Lysine-Trimethylation", "K", "+42.047", false, new ChemicalForm(3, 6, 0, 0, 0));
		Modification TRIMETHYL_R = new Modification("Arginine-Trimethylation", "R", "+42.047", false, new ChemicalForm(3, 6, 0, 0, 0));
		Modification TETRAMETHYL_K = new Modification("Lysine-Tetramethylation", "K", "+56.063", false, new ChemicalForm(4, 8, 0, 0, 0));
		Modification TETRAMETHYL_R = new Modification("Arginine-Tetramethylation", "R", "+56.063", false, new ChemicalForm(4, 8, 0, 0, 0));

		Modification ACETHYL_K = new Modification("Lysine-Acethylation", "K", "+42.011", false, new ChemicalForm(2, 2, 1, 0, 0));
		Modification ACETHYL_NTERM = new Modification("Nterm-Acethylation", ANY_NTERM, "+42.011", false, new ChemicalForm(2, 2, 1, 0, 0));

		Modification AMMONIALOSS_N = new Modification("Asparagine-Ammonialoss", "N", "-17.027", false, new ChemicalForm(0, -3, 0, -1, 0));
		Modification DEHYDRATION_S = new Modification("Serine-Dehydration", "S", "-18.011", false, new ChemicalForm(0, -2, -1, 0, 0));
		Modification DEHYDRATION_T = new Modification("Threonine-Dehydration", "T", "-18.011", false, new ChemicalForm(0, -2, -1, 0, 0));
		Modification DEHYDRATION_D = new Modification("AsparticAcid-Dehydration", "D", "-18.011", false, new ChemicalForm(0, -2, -1, 0, 0));
		Modification DEHYDRATION_C = new Modification("Cysteine-Dehydration", ANY_NTERM, "-18.011", false, new ChemicalForm(0, -2, -1, 0, 0));

		Modification DEAMIDATION_N = new Modification("Asparagine-Deamidation", "N", "+0.984", false, new ChemicalForm(0, -1, 1, -1, 0));
		Modification DEAMIDATION_Q = new Modification("Glutamine-Deamidation", "Q", "+0.984", false, new ChemicalForm(0, -1, 1, -1, 0));
		Modification CARBON_NTERM = new Modification("Nterm-Carbon", ANY_NTERM, "+12.000", false, new ChemicalForm(1, 0, 0, 0, 0));

		Modification PHOSPHORYLATION_T = new Modification("Threonine-Phosphorylation", "T", "+79.966", false, new ChemicalForm(0, 1, 3, 0, 0, 1));
		Modification PHOSPHORYLATION_Y = new Modification("Tyroshine-Phosphorylation", "Y", "+79.966", false, new ChemicalForm(0, 1, 3, 0, 0, 1));
		Modification PHOSPHORYLATION_S = new Modification("Serine-Phosphorylation", "S", "+79.966", false, new ChemicalForm(0, 1, 3, 0, 0, 1));

		Modification FORMYL_NTERM = new Modification("Nterm-Formyl", ANY_NTERM, "+27.995", false, new ChemicalForm(1, 0, 1, 0, 0));
		Modification FORMYL_T = new Modification("Threonine-Formyl", "T", "+27.995", false, new ChemicalForm(1, 0, 1, 0, 0));
		Modification FORMYL_S = new Modification("Serine-Formyl", "S", "+27.995", false, new ChemicalForm(1, 0, 1, 0, 0));

		Modification AEBS_K = new Modification("Lysine-AEBS", "K", "+183.035", false, new ChemicalForm(8, 9, 2, 1, 1));
		Modification AEBS_Y = new Modification("Tyroshine-AEBS", "Y", "+183.035", false, new ChemicalForm(8, 9, 2, 1, 1));

		Modification DETHIOEMTHYL_M = new Modification("Methionine-Dethiomethyl", "M", "-48.003", false, new ChemicalForm(-1, -4, 0, 0, -1));
		Modification DEHYDROALANINE_C = new Modification("Cysteine-Dehydroalanine", "C", "-33.988", false, new ChemicalForm(0, -2, 0, 0, -1));

		Modification SULFURDIOXIDE_C = new Modification("Cysteine-SulfurDioxide", "C", "+63.962", false, new ChemicalForm(0, 0, 2, 0, 1));

		Modification GLYGLY_C = new Modification("Cysteine-GlyGly", "C", "+114.043", false, new ChemicalForm(4, 6, 2, 2, 0));
		Modification GLYGLY_K = new Modification("Lysine-GlyGly", "K", "+114.043", false, new ChemicalForm(4, 6, 2, 2, 0));
		Modification GLYGLY_T = new Modification("Threonine-GlyGly", "T", "+114.043", false, new ChemicalForm(4, 6, 2, 2, 0));
		Modification GLYGLY_S = new Modification("Serine-GlyGly", "S", "+114.043", false, new ChemicalForm(4, 6, 2, 2, 0));

		Modification CYSSER_C = new Modification("Cysteine-CysSer", "C", "-15.977", false, new ChemicalForm(0, 0, 1, 0, -1));
		Modification NITROSYL_C = new Modification("Cysteine-Nitrosyl", "C", "+28.990", false, new ChemicalForm(0, -1, 1, 1, 0));


		Modification TRPKYNURENIN_W = new Modification("Tryptophan-TrpKynurenin", "W", "+3.995", false, new ChemicalForm(-1, 0, 1, 0, 0));

		chemicalForm = new ChemicalForm(0, -1, 0, 0, 0);
		chemicalForm.setI(1);
		Modification IODO_Y = new Modification("Tyroshine-Iodo", "Y", "+125.897", false, chemicalForm);

		addModification(TMT_NTERM);
		addModification(TMT_K);
		addModification(CARBAM_NTERM);
		addModification(CARBAM_C);
		addModification(CARBAM_K);
		addModification(CARBAM_H);
		addModification(CARBAMDTT_C);
		addModification(CARBAMYL_K);
		addModification(CARBAMYL_NTERM);
		addModification(CATION_NA_E);
		addModification(CATION_NA_D);
		addModification(CATION_FE_II_D);
		addModification(CATION_FE_II_E);
		addModification(CATION_K_E);
		addModification(CATION_K_D);
		addModification(OXIDATION_M);
		addModification(OXIDATION_P);
		addModification(OXIDATION_W);
		addModification(OXIDATION_C);
		addModification(DIOXIDATION_M);
		addModification(DIOXIDATION_C);
		addModification(DIOXIDATION_W);
		addModification(TRIOXIDATION_C);
		addModification(METHYL_K);
		addModification(METHYL_R);
		addModification(METHYL_E);
		addModification(METHYL_D);
		addModification(METHYL_H);
		addModification(DIMETHYL_K);
		addModification(DIMETHYL_R);
		addModification(TRIMETHYL_R);
		addModification(TRIMETHYL_K);
		addModification(TETRAMETHYL_K);
		addModification(TETRAMETHYL_R);
		addModification(ACETHYL_K);
		addModification(ACETHYL_NTERM);
		addModification(AMMONIALOSS_N);
		addModification(DEHYDRATION_S);
		addModification(DEHYDRATION_D);
		addModification(DEHYDRATION_T);
		addModification(DEAMIDATION_N);
		addModification(DEAMIDATION_Q);
		addModification(DEHYDRATION_C);
		addModification(PHOSPHORYLATION_T);
		addModification(PHOSPHORYLATION_Y);
		addModification(PHOSPHORYLATION_S);
		addModification(CARBON_NTERM);
		addModification(FORMYL_NTERM);
		addModification(FORMYL_T);
		addModification(FORMYL_S);
		addModification(AEBS_K);
		addModification(AEBS_Y);
		addModification(DETHIOEMTHYL_M);
		addModification(DEHYDROALANINE_C);
		addModification(SULFURDIOXIDE_C);
		addModification(GLYGLY_S);
		addModification(GLYGLY_T);
		addModification(GLYGLY_K);
		addModification(GLYGLY_C);
		addModification(CYSSER_C);
		addModification(NITROSYL_C);
		addModification(TRPKYNURENIN_W);
		addModification(IODO_Y);
	}

	public static void addModification (Modification mod) {
		if(preDefinedModificationList == null) {
			preDefinedModificationList = new Hashtable<>();
			massByModName = new Hashtable<>();
			initModification();
		}

		ArrayList<Modification> mods = preDefinedModificationList.get(mod.modSite+mod.modMass);
		if(mods == null) {
			mods = new ArrayList<>();
		}
		mods.add(mod);

		preDefinedModificationList.put(mod.modSite+mod.modMass, mods);
		massByModName.put(mod.modName, mod.modMass);
	}

	public ChemicalForm getChemicalForm () {
		return this.chemicalForm;
	}

	public double getModMass () {
		if(getChemicalForm() != null) {
			return getChemicalForm().getMass();
		} else {
			return Double.parseDouble(this.modMass);
		}
	}

	public String getModName () {
		return this.modName;
	}

	public String getModSite () {
		return this.modSite;
	}

	public boolean isFixed () {
		return this.isFixed;
	}

	public static ArrayList<Modification> getListOfModifcation () {
		ArrayList<Modification> modifications = new ArrayList<>();
		Iterator<String> keys = (Iterator<String>)Modification.preDefinedModificationList.keys();
		while(keys.hasNext()) {
			ArrayList<Modification> modification = Modification.preDefinedModificationList.get(keys.next());
			modifications.addAll(modification);
		}
		return modifications;
	}

}