# HYCOM demo

## HYCOM version

`v2.3.01`

## How-to instructions

1. Compile HYCOM and obtain its `compile_commands.json` file:

```bash
$ cd src && bear -- make ARCH=codee TYPE=demo
```

2. Use `codee screening` to obtain a high-level overview of the project:

```bash
$ codee screening --include-categories all --config compile_commands.json
```

3. Use `codee checks --verbose` to further investigate the checkers
identified by Codee:

```bash
$ codee checks --verbose --include-categories all --config compile_commands.json
```

4. Use `codee checks --verbose --level 1` to focus on the most severe
checkers identified by Codee:

```bash
$ codee checks --verbose --level 1 --include-categories all --config compile_commands.json
```
