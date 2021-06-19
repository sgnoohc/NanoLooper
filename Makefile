SUBDIRS := NanoLooper
LIBRARIES := NanoTools/NanoCORE rooutil

all: $(LIBRARIES) $(SUBDIRS)

NanoCORE:
	$(MAKE) -C NanoTools/NanoCORE

rooutil:
	$(MAKE) -C rooutil

$(SUBDIRS): NanoCORE rooutil
	$(MAKE) -C $@

.PHONY: all $(LIBRARIES) $(SUBDIRS)

clean:
	cd NanoLooper/ && make clean;

cleanall:
	cd rooutil/ && make clean;
	cd NanoTools/NanoCORE/ && make clean;
	cd NanoLooper/ && make clean;
