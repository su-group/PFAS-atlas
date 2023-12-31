import torch
from rxnfp.models import SmilesLanguageModelingModel
import os


os.environ['CUDA_LAUNCH_BLOCKING'] = '1'

config = {
    "architectures": [
        "BertForMaskedLM"
    ],
    "attention_probs_dropout_prob": 0.1,
    "hidden_act": "gelu",
    "hidden_dropout_prob": 0.1,
    "hidden_size": 256,
    "initializer_range": 0.02,
    "intermediate_size": 512,
    "layer_norm_eps": 1e-12,
    "max_position_embeddings": 512,
    "model_type": "bert",
    "num_attention_heads": 4,
    "num_hidden_layers": 12,
    "pad_token_id": 0,
    "type_vocab_size": 2,

}

vocab_path = 'Bert_model/vocab.txt'

args = {'config': config,
        'vocab_path': vocab_path,
        'train_batch_size': 16,
        'manual_seed': 42,
        'fp16': False,
        "num_train_epochs": 10,
        'max_seq_length': 300,
        'evaluate_during_training': True,
        'overwrite_output_dir': True,
        'output_dir': 'out/mlm_QM100w',
        'learning_rate': 1e-4,
        # 'wandb_project': "pretrain_exp_0708"
        }

model = SmilesLanguageModelingModel(model_type='bert', model_name=None, args=args, use_cuda=torch.cuda.is_available())

train_file = 'Bert_model/step3_OECD_Class_0812.csv'
eval_file = 'Bert_model/step3_OECD_Class_0812.csv'

model.train_model(train_file=train_file, eval_file=eval_file)
