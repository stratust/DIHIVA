script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
pipeline_dir=$(echo $script_dir | perl -plane 's/^(.*)\/.*/$1/g')
dir=$1

ln -s ${pipeline_dir}/Snakefile ${dir}/local.Snakefile
ln -s ${pipeline_dir}/config $dir
ln -s ${pipeline_dir}/config.yml $dir
ln -s ${pipeline_dir}/database $dir
ln -s ${pipeline_dir}/scripts $dir
ln -s ${pipeline_dir}/subworkflows $dir
ln -s ${pipeline_dir}/shell_scripts/run_local_no_cluster.sh $dir
