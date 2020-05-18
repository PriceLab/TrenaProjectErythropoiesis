ID = 18may2020
IMAGE_NAME = srm_rna_combined$(ID)
TESTPORT=60055

build:
	docker build -t $(IMAGE_NAME) -f Dockerfile .

bash:
	docker run -ti --rm -p $(TESTPORT):3838 $(IMAGE_NAME) bash

shinyRun:
	docker run -ti --rm -p $(TESTPORT):3838 rocker/shiny bash

dockerRun:
	docker run -p $(TESTPORT):3838 $(IMAGE_NAME)

save:
	time docker save $(IMAGE_NAME) >  ~/tmp/docker_$(IMAGE_NAME).image
	ls -lat ~/tmp/docker_$(IMAGE_NAME).image

scp:
	time scp ~/tmp/docker_$(IMAGE_NAME).image pshannon@trena:/local/tmp/
	ssh pshannon@trena ls -l /local/tmp/docker_$(IMAGE_NAME).image
#	ssh pshannon@trena docker load < /local/tmp/docker_$(IMAGE_NAME).image


