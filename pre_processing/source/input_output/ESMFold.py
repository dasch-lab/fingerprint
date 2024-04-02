import torch
from transformers import AutoTokenizer, EsmForProteinFolding
from transformers.models.esm.openfold_utils.protein import to_pdb, Protein as OFProtein
from transformers.models.esm.openfold_utils.feats import atom14_to_atom37

# Create the model
is_cuda = torch.cuda.is_available()
# model = esm.pretrained.esmfold_v1()
# model = model.eval().cuda()
from Bio.PDB import PDBParser
from Bio import SeqIO

esm_model = None
tokenizer = None

def init_model():

  global esm_model
  global tokenizer
  
  if esm_model is not None:
    return
  
  tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
  esm_model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
  if is_cuda:
    esm_model = esm_model.cuda()

  # enable TensorFloat32 computation for a general speedup
  torch.backends.cuda.matmul.allow_tf32 = True

  # Uncomment this line if your GPU memory is 16GB or less, or if you're folding longer (over 600 or so) sequences
  esm_model.trunk.set_chunk_size(64)


class ESMModel():

  def generate_model(self, model_path, data):

    # https://github.com/huggingface/notebooks/blob/main/examples/protein_folding.ipynb

    # Make sure the model is initalized
    init_model()
    
    # If you're using a GPU, you'll need to move the tokenized data to the GPU now.
    tokenized_input = tokenizer(data, return_tensors="pt", add_special_tokens=False)['input_ids']
    if is_cuda:
      tokenized_input = tokenized_input.cuda()

    # Predict structure
    with torch.no_grad():
      output = esm_model(tokenized_input)

    # Extract the model structure
    pdb = ESMModel.convert_outputs_to_pdb(output)

    # Store result to file
    with open(model_path, "w") as handle:
      handle.write("".join(pdb))

    return True
  
  @staticmethod
  def convert_outputs_to_pdb(outputs):
    final_atom_positions = atom14_to_atom37(outputs["positions"][-1], outputs)
    outputs = {k: v.to("cpu").numpy() for k, v in outputs.items()}
    final_atom_positions = final_atom_positions.cpu().numpy()
    final_atom_mask = outputs["atom37_atom_exists"]
    pdbs = []
    for i in range(outputs["aatype"].shape[0]):
      aa = outputs["aatype"][i]
      pred_pos = final_atom_positions[i]
      mask = final_atom_mask[i]
      resid = outputs["residue_index"][i] + 1
      pred = OFProtein(
          aatype=aa,
          atom_positions=pred_pos,
          atom_mask=mask,
          residue_index=resid,
          b_factors=outputs["plddt"][i],
          chain_index=outputs["chain_index"][i] if "chain_index" in outputs else None,
      )
      pdbs.append(to_pdb(pred))       
      
    return pdbs
  
  @staticmethod
  def get_b_factor():
     pass

"""
def parse_pdb_b_factors(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('pdb_structure', pdb_file)
    
    b_factors = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    b_factors.append(atom.get_bfactor())
    
    return b_factors


if __name__ == "__main__":
  #parser = PDBParser(QUIET=True)
  #structure = parser.get_structure('structure', "/disk1/fingerprint/tmp/1AKJ_AB.pdb")
  sequence = []
  for record in SeqIO.parse("/disk1/fingerprint/tmp/4FQI_HL.pdb", "pdb-atom"):
    print(record.id[-1])
    sequence.append(str(record.seq))

  model = ESMModel()
  model.generate_model("/disk1/fingerprint/provaESMFold/output.pdb",sequence[0])

  pdb_data = parse_pdb_b_factors("/disk1/fingerprint/provaESMFold/output.pdb")
  print(pdb_data)
"""