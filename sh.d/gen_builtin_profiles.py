import base64
import os

sname_file_mapping = {
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


if __name__ == "__main__":
    with (
        open(os.path.join("src", "art", "builtin_profiles.cc"), "w") as w,
        open(os.path.join("src", "art", "builtin_profiles.hh"), "w") as wh,
    ):
        wh.write("#pragma once\n")
        wh.write("namespace labw::art_modern {\n")
        w.write("namespace labw::art_modern {\n")
        wh.write(f'char NULL_PROFILE[] = "\\0";\n')
        snames = []
        snames_constructed = []
        for sname, files in sname_file_mapping.items():
            snames.append(sname)
            for i, file in enumerate(files):
                with open(os.path.join("data", "Illumina_profiles", file), "rb") as r:
                    data = str(base64.b64encode(r.read()), encoding="US-ASCII")

                w.write(f"char {sname}_{i}[] = \\")
                while data:
                    w.write('\n"')
                    w.write(data[:80])
                    data = data[80:]
                    w.write('"')
                w.write(";\n")

                wh.write(f"extern char {sname}_{i}[];\n")
            snames_constructed.append("{" + f"{sname}_0, " + ("NULL_PROFILE" if i == 0 else f"{sname}_1") + "}")
        wh.write(f"const int N_BUILTIN_PROFILE = {len(snames_constructed)};\n")
        wh.write("char* ENCODED_BUILTIN_PROFILES[N_BUILTIN_PROFILE][2] = {\n")
        for sname_constructed in snames_constructed:
            wh.write(f"    {sname_constructed},\n")
        wh.write("};\n")

        wh.write("const char* const BUILTIN_PROFILE_NAMES[N_BUILTIN_PROFILE] = {\n")
        for sname in snames:
            wh.write(f'    "{sname}",\n')
        wh.write("};\n")
        wh.write("} // namespace labw::art_modern\n")
        w.write("} // namespace labw::art_modern\n")
