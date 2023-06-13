name=$1
version=$2

aws ecr-public get-login-password --region us-east-1 | docker login --username AWS --password-stdin public.ecr.aws/csgenetics

docker build -t $name -t $name:$version --platform linux/x86_64 $name/.
docker tag $name:latest public.ecr.aws/n1x9b6b5/$name:latest
docker tag $name:$version public.ecr.aws/n1x9b6b5/$name:$version
docker push public.ecr.aws/n1x9b6b5/$name:latest
docker push public.ecr.aws/n1x9b6b5/$name:$version
