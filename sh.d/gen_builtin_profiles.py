import datetime
import gzip
import os
import sys

SNAME_FILE_MAPPING = {
    # TODO: Separate this data to some INI/JSON file
    "GA1_36bp": ("Emp36R1.txt", "Emp36R2.txt"),
    "GA1_44bp": ("Emp44R1.txt", "Emp44R2.txt"),
    "GA2_50bp": ("Emp50R1.txt", "Emp50R2.txt"),
    "GA2_75bp": ("Emp75R1.txt", "Emp75R2.txt"),
    "HiSeq1000_100bp": ("Emp100R1.txt", "Emp100R2.txt"),
    "MiSeq_250bp": ("EmpMiSeq250R1.txt", "EmpMiSeq250R2.txt"),
    "GA1Recalibrated_36bp": ("EmpR36R1.txt", "EmpR36R2.txt"),
    "GA1Recalibrated_44bp": ("EmpR44R1.txt", "EmpR44R2.txt"),
    "GA2Recalibrated_50bp": ("EmpR50R1.txt", "EmpR50R2.txt"),
    "GA2Recalibrated_75bp": ("EmpR75R1.txt", "EmpR75R2.txt"),
    "HiSeq2000_100bp": ("HiSeq2kL100R1.txt", "HiSeq2kL100R2.txt"),
    "HiSeq2500_125bp": ("HiSeq2500L125R1.txt", "HiSeq2500L125R2.txt"),
    "HiSeq2500_150bp": ("HiSeq2500L150R1.txt", "HiSeq2500L150R2.txt"),
    "HiSeq2500Filtered_150bp": ("HiSeq2500L150R1filter.txt", "HiSeq2500L150R2filter.txt"),
    "HiSeqX_PCR_Free_150bp": ("HiSeqXPCRfreeL150R1.txt", "HiSeqXPCRfreeL150R2.txt"),
    "HiSeqX_TruSeq_150bp": ("HiSeqXtruSeqL150R1.txt", "HiSeqXtruSeqL150R2.txt"),
    "MiniSeq_TruSeq_50bp": ("MiniSeqTruSeqL50.txt",),
    "MiSeq_v3_250bp": ("MiSeqv3L250R1.txt", "MiSeqv3L250R2.txt"),
    "NextSeq500_v2_75bp": ("NextSeq500v2L75R1.txt", "NextSeq500v2L75R2.txt"),
}

DTYPE = "unsigned char"
AUTOGEN_HEADER = (
    f"// This file is auto-generated "
    f"by {os.path.basename(__file__)} "
    f"at {datetime.datetime.now().isoformat()}\n"
)
NAMESPACE_HEADER = "namespace labw::art_modern {\n"
NAMESPACE_FOOTER = "}\n // namespace labw::art_modern\n"

if __name__ == "__main__":
    os.makedirs(sys.argv[1], exist_ok=True)
    with open(os.path.join(sys.argv[1], "builtin_profiles.cc"), "w", encoding="US-ASCII") as wc:
        with open(os.path.join(sys.argv[1], "builtin_profiles.hh"), "w", encoding="US-ASCII") as wh:
            wh.write(AUTOGEN_HEADER)
            wc.write(AUTOGEN_HEADER)

            wh.write("#pragma once\n")
            wh.write("#include <cstddef>\n")
            wh.write(NAMESPACE_HEADER)
            wc.write(NAMESPACE_HEADER)

            wh.write(f"extern {DTYPE} NULL_PROFILE[1];\n")
            wc.write(DTYPE + " NULL_PROFILE[1] = {0};\n")

            snames = []
            snames_constructed = []
            slengths_constructed = []
            for sname, files in SNAME_FILE_MAPPING.items():
                slengths = []
                snames.append(sname)
                for i, file in enumerate(files):
                    with open(os.path.join("data", "Illumina_profiles", file), "rb") as r:
                        data = gzip.compress(r.read(), 9)

                    with open(os.path.join(sys.argv[1], f"{sname}_{i}.cc"), "w", encoding="US-ASCII") as sw:
                        sw.write(AUTOGEN_HEADER)
                        sw.write(NAMESPACE_HEADER)

                        sw.write(f"{DTYPE} {sname}_{i}[{len(data)}]" + " = {\n")
                        cw = 0
                        for c in data[:-1]:
                            sw.write(hex(c))
                            sw.write(", ")
                            cw += 1
                            if cw % 10 == 0:
                                sw.write("\n")
                        sw.write(hex(data[-1]))
                        sw.write("\n")
                        sw.write("};\n")
                        sw.write(NAMESPACE_FOOTER)
                    slengths.append(len(data))

                    wh.write(f"extern {DTYPE} {sname}_{i}[{len(data)}];\n")
                snames_constructed.append("{" + f"{sname}_0, " + ("NULL_PROFILE" if i == 0 else f"{sname}_1") + "}")
                slengths_constructed.append("{" + f"{slengths[0]}, " + ("0" if i == 0 else f"{slengths[1]}") + "}")

            wh.write(f"const int N_BUILTIN_PROFILE = {len(snames_constructed)};\n")
            wh.write(DTYPE + "* ENCODED_BUILTIN_PROFILES[N_BUILTIN_PROFILE][2] = {\n")
            for sname_constructed in snames_constructed:
                wh.write(f"    {sname_constructed},\n")
            wh.write("};\n")

            wh.write("const char* const BUILTIN_PROFILE_NAMES[N_BUILTIN_PROFILE] = {\n")
            for sname in snames:
                wh.write(f'    "{sname}",\n')
            wh.write("};\n")

            wh.write("const std::size_t BUILTIN_PROFILE_LENGTHS[N_BUILTIN_PROFILE][2] = {\n")
            for slength_constructed in slengths_constructed:
                wh.write("    " + slength_constructed + ",\n")
            wh.write("};\n")

            wh.write(NAMESPACE_FOOTER)
            wc.write(NAMESPACE_FOOTER)
