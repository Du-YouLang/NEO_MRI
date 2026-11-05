"""HD-BET library adapter: CPU4, noTTA, one preprocessing/export worker."""
def main(argv=None):
    import argparse,os
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',required=True);p.add_argument('--output',required=True);p.add_argument('--model-dir',required=True);a=p.parse_args(argv)
    os.environ['CUDA_VISIBLE_DEVICES']='-1';os.environ['nnUNet_compile']='F'
    if hasattr(os,'sched_setaffinity'):os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[:4])
    import HD_BET.paths as paths
    paths.folder_with_parameter_files=a.model_dir
    import torch
    from HD_BET.hd_bet_prediction import get_hdbet_predictor
    predictor=get_hdbet_predictor(use_tta=False,device=torch.device('cpu'))
    # Upstream resets CPU threads to os.cpu_count(); enforce the reviewed budget afterwards.
    torch.set_num_threads(4)
    predictor.predict_from_files([[a.input]],[a.output],save_probabilities=False,overwrite=False,
        num_processes_preprocessing=1,num_processes_segmentation_export=1)
    return 0
if __name__=='__main__':raise SystemExit(main())
