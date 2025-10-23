---
permalink: /research/
title: "Community Research"
excerpt: "Contributing to iFEM and show your research"
sidebar:
    nav: research
---

## Contributing

We welcome your input and feedback on the development of iFEM. You can

- Reporting a bug (issue)
- Submitting a fix (pull request)
- Proposing new features (pull request or issue)
- Show your research using iFEM (pull request, see below)

## Community Research

To show your research using iFEM, you can submit a pull request with the following files:
- In `/research/` folder, make a new folder with the name of your research, e.g. `/research/my_research`. Put a minimal working example in this folder that is runnable under the current iFEM version.
- In `/docs/_research` folder, write a `YOURNAME_RESEARCHNAME.md` file with the following content
    * Introduction of you and your research.
    * The core code modified or created under iFEM framework in this research.
    * (Optional) Bibliography of your research.
- (Optional) In `/docs/_data/author.yml`, add your name to the `authors` list following the template.
- Submit a pull request at [iFEM GitHub](https://github.com/lyc102/ifem).