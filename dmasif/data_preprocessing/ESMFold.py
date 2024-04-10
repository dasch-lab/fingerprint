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
"""
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

  def generate_model(self, data, pdb_write = True, model_path = None):

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
    
    if pdb_write:
      # Extract the model structure
      pdb = ESMModel.convert_outputs_to_pdb(output)

      # Store result to file
      with open(model_path, "w") as handle:
        handle.write("".join(pdb))
      
      return True
    else:
       outputs = {k: v.to("cpu").numpy() for k, v in outputs.items()}
       return output["plddt"]
  
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
  
  def get_plddt(output):
    outputs = {k: v.to("cpu").numpy() for k, v in outputs.items()}
    return output["plddt"]
     

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

def dictionary_pdb_bfactors(structure):
   b_factor_dict = []
   for model in structure:
        for chain in model:
            chain =  chain.get_id()
            for residue in chain:
                residue_id = residue.get_id()[1]
                residue_name = residue.get_resname()
                for atom in residue:
                    atom_name =  atom.get_name() 
                    key = chain + '_' + residue_id + '_' + residue_name + '_' + atom_name
                    b_factor_dict[key] = atom.get_factor()
   return b_factor_dict

if __name__ == "__main__":
  #parser = PDBParser(QUIET=True)
  #structure = parser.get_structure('structure', "/disk1/fingerprint/tmp/1AKJ_AB.pdb")
  sequence = []
  for record in SeqIO.parse("/disk1/fingerprint/tmp/4FQI_HL.pdb", "pdb-atom"):
    print(record.id[-1])
    sequence.append(str(record.seq))

  model = ESMModel()
  output = model.generate_model(data = sequence[0], pdb_write = True, model_path = "/disk1/fingerprint/provaESMFold/output.pdb")

  pdb_data = parse_pdb_b_factors("/disk1/fingerprint/provaESMFold/4FQI_HL_ESMFold.pdb")
  print(pdb_data)"""

import torch
from transformers import AutoTokenizer, EsmForProteinFolding
from transformers.models.esm.openfold_utils.protein import to_pdb, Protein as OFProtein
from transformers.models.esm.openfold_utils.feats import atom14_to_atom37

# Create the model
is_cuda = False  # Set CUDA availability to False

# Initialize the model and tokenizer
def init_model():
    global esm_model
    global tokenizer
    
    if esm_model is not None:
        return
    
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    esm_model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
    
    # No need to move the model to CUDA
    # No need to set CUDA-specific configurations

    # Uncomment this line if your GPU memory is 16GB or less, or if you're folding longer (over 600 or so) sequences
    esm_model.trunk.set_chunk_size(64)

class ESMModel():
    def generate_model(self, chain, data, pdb_write=True, model_path=None):
        init_model()
        plddt = []
        
        # Tokenized input is on CPU by default
        tokenized_input = tokenizer(data, return_tensors="pt", add_special_tokens=False)['input_ids']
    
        # Predict structure
        with torch.no_grad():
            output = esm_model(tokenized_input)
        
        if pdb_write:
            # Extract the model structure
            pdb = ESMModel.convert_outputs_to_pdb(output, chain)
            
            # Store result to file
            with open(model_path, "w") as handle:
              pdb_with_modified_chain = ESMModel.modify_chain_id(pdb, chain)
              handle.write(pdb_with_modified_chain)

              #handle.write("".join(pdb))
              return True
        else:
            output = {k: v.to("cpu").numpy() for k, v in output.items()}
            plddt.append(output["plddt"])
            return plddt
            
        
    @staticmethod
    def convert_outputs_to_pdb(outputs, chain):
        # Convert tensors to NumPy arrays on CPU
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
    
    def modify_chain_id(pdb, new_chain_id):
      lines = pdb[0].split('\n')
      modified_lines = []
      for line in lines:
          if line.startswith('ATOM') or line.startswith('HETATM'):
              line = line[:21] + new_chain_id + line[22:]
          modified_lines.append(line)
      return '\n'.join(modified_lines) + '\n'

    def get_plddt(output):
        outputs = {k: v.to("cpu").numpy() for k, v in outputs.items()}
        print(output)
        return output["plddt"]

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

def dictionary_pdb_bfactors(structure, b_factor_dict = None):
   if b_factor_dict is None:
       b_factor_dict = {}
   
   for atom in structure.get_atoms():
       atom_name = atom.get_name()
       residue_name = atom.get_parent().get_resname()
       chain = atom.get_parent().get_parent().id
       residue_id = atom.get_parent().get_id()[1]
       b_factor = atom.get_bfactor()
       key = chain + '_' + str(residue_id) + '_' + residue_name + '_' + atom_name
       b_factor_dict[key] = atom.get_bfactor()
   return b_factor_dict

if __name__ == "__main__":
    model = ESMModel()
    sequence = {}
    for record in SeqIO.parse("/disk1/fingerprint/tmp/4FQI_AB.pdb", "pdb-atom"):
        #print(record.id[-1])
        chain = record.id[-1]
        sequence[record.id[-1]] = str(record.seq)
        out_file = '4FQI_'+ chain + '_ESMFold.pdb'
        path = '/'.join(["/disk1/fingerprint/provaESMFold", out_file])
        output = model.generate_model(chain = chain, data=sequence[chain], pdb_write=True, model_path=path)
