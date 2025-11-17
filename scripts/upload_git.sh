=================================
#Upload project to github
==============================
#1.install DVC

conda activate bioinfo

if command -v mamba >/dev/null 2>&1; then
  mamba install -y -c conda-forge dvc dvc-ssh dvc-s3
else
  conda install -y -c conda-forge dvc dvc-ssh dvc-s3
fi

hash -r
dvc --version
 

 #2.INITIALIZE GIT AND DVC

 cd /Users/mm/Documents/project/metagenome
# init git if needed
[ -d .git ] || git init
# init dvc if needed
[ -d .dvc ] || dvc init && git add .dvc .dvcignore && git commit -m "init dvc"


#3.track large data with DVC

#FIND PATH
cd /Users/mm/Documents/project/metagenome
ls -la
# search for the MAG folder(s)
find . -maxdepth 3 -type d -name "mag_*" -print
# or search for any "bins" dirs
find . -type d -name "bins" -print

#ADD TO DVC
dvc add 1_quality_control/trimmed_reads || true 
dvc add ./3_MAG_gene_cluster/mag_SRR24442557/bins || true 
dvc add ./3_MAG_gene_cluster/mag_SRR24442557/assembly || true
dvc add ./0_rawdata || true 

#raw data is still in git. letes safe migrate it to dvc

# stop Git from tracking the raw data (keep files on disk)
git rm -r --cached 0_rawdata
git commit -m "stop tracking 0_rawdata in git (move to DVC)"

# add to DVC and commit pointer
dvc add 0_rawdata
git add 0_rawdata.dvc .gitignore
git commit -m "dvc: track 0_rawdata"
dvc add ./0_rawdata || true 

brew install awscli
aws --version
aws configure
aws s3 ls

#create s3 bucket
aws s3 mb s3://my-metagenome-dvc --region ap-southeast-1
#result: make_bucket: my-metagenome-dvc
aws s3 ls #comfirmed: my-metagenome-dvc

#tell dvc to use this bucket
dvc remote add -d storage s3://my-metagenome-dvc

#it seems i have old remote. lets delete the old remote
dvc remote remove storage
#add new remote
dvc remote add -d storage s3://my-metagenome-dvc
dvc push

#the process keep restarting
dvc status -c
dvc add 0_rawdata
git add 0_rawdata.dvc .gitignore
dvc status -c
dvc list . --dvc-only
dvc list 0_rawdata #weird result

dvc list . --dvc-only


#check if the files is uploaded in DVC
dvc status -c
dvc list . --dvc-only
dvc list s3://my-metagenome-dvc --recursive
dvc push -v


# Check that everything is synced
dvc status -c

# Then clean local cache safely
dvc gc -w


###lets upload to github
git remote -v
cat .gitignore
# Ignore DVC cache
.dvc/cache/
ls -lh .dvc/cache/
git check-ignore -v .dvc/cache/

git add .dvc/config *.dvc .gitignore
git commit -m "Add DVC config and tracked data files"
git push origin main

git pull origin main --no-rebase

git pull origin main --allow-unrelated-histories #failed, 16 gb in total


# Undo big push attempt (doesn’t delete work)
git reset --soft HEAD~1

# Re-add only small necessary files
git add .dvc/config *.dvc .gitignore *.sh README.md

# Commit again
git commit -m "Track data with DVC (data stored in S3 remote)"

# Push safely
git push origin main #fail


# 1️⃣ Fetch the latest remote commits
git fetch origin main

# 2️⃣ Merge them safely into your local branch
git pull origin main --allow-unrelated-histories --no-rebase


git add .

rm -rf 2_taxonomy_functional/bracken/.git
git add 2_taxonomy_functional/bracken
git commit -m "Remove embedded Bracken git repo and finalize merge"
git push origin main
#notworking

git lfs ls-files
du -h --max-depth=2 | sort -hr | head -20








# Ignore DVC temporary/lock files
*.dvc.lock








