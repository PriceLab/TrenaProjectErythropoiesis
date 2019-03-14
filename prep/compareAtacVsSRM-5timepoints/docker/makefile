ID = 14mar2019
IMAGE_NAME = srm_atac_fli1_klf1_$(ID)
TESTPORT=60055

build:
	docker build -t $(IMAGE_NAME) -f Dockerfile .

bash:
	docker run -ti --rm -p $(TESTPORT):3838 $(IMAGE_NAME) bash

dockerRun:
	docker run -p $(TESTPORT):3838 $(IMAGE_NAME)

save:
	time docker save $(IMAGE_NAME) >  ~/tmp/docker_$(IMAGE_NAME).image
	ls -lat ~/tmp/docker_$(IMAGE_NAME).image

scp:
	time scp ~/tmp/docker_$(IMAGE_NAME).image pshannon@trena:/local/tmp/
	ssh pshannon@trena ls -l /local/tmp/docker_$(IMAGE_NAME).image
#	ssh pshannon@trena docker load < /local/tmp/docker_$(IMAGE_NAME).image


