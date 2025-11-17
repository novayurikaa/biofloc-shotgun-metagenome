

cd /Users/mm/Documents/project/metagenome
git init
du -h -d 3 .

#add folder structure first
#1. raw data
cd /Users/mm/Documents/project/metagenome/0_rawdata
touch .gitkeep
git add -f .gitkeep. #we use force since rawdata is in gitignore
git commit -m "Add empty 0_rawdata folder for folder structure"
git push origin main


