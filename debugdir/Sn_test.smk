rule all:
    input: 
        "debugcmd.sh", "test1.txt", "test2.txt"
    



rule test1:
    input: "debugcmd.sh"
    output: 
        one = "test1.txt"
    shell : "STAR --help > {output.one}"

rule test2:
    input: rules.test1.output.one
    output: 
        two = "test2.txt"
    shell: "echo END > {output.two}"
