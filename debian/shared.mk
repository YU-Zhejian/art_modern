
ifeq ($(wildcard /usr/include/concurrentqueue/moodycamel/concurrentqueue.h),/usr/include/concurrentqueue/moodycamel/concurrentqueue.h)
    USE_CONCURRENT_QUEUE_PATH := /usr/include/concurrentqueue/moodycamel/
else ifeq ($(wildcard /usr/include/concurrentqueue/concurrentqueue.h),/usr/include/concurrentqueue/concurrentqueue.h)
    USE_CONCURRENT_QUEUE_PATH := /usr/include/concurrentqueue/
else
    $(error concurrentqueue.h not found in expected locations.)
endif
$(info Using moodycamel path: $(USE_CONCURRENT_QUEUE_PATH))

%.1: %.1.md
	# TODO: Supress: art-modern: wrong-manual-section 1 != 7 [usr/share/man/man1/art_modern.1.gz:2]
	# Tried: sed -E 's;^\.TH "([^"]+)" "([^"]+)" "([^"]+)";.TH \1 \2 \3;'
	lowdown -s -Tman -o /dev/stdout $< > $@

%:
	dh $@ --buildsystem=cmake

override_dh_makeshlibs:
# We do not install libraries in standard locations, so we skip this step
