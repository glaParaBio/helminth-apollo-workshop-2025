<!-- vim-markdown-toc GFM -->

* [Setting up](#setting-up)
    * [Getting the data](#getting-the-data)
    * [Running Apollo](#running-apollo)

<!-- vim-markdown-toc -->

# Setting up

To run this demo, you'll need to have
[Docker](//docs.docker.com/engine/install/) (minimum version 27) and
[Docker Compose](//docs.docker.com/compose/install/) (minimum version 2.17.0)
installed, as well as [`jq`](https://jqlang.github.io/jq/).

We also recommend a Node.js package manager (such as `npm` or `yarn`) to make
things a easier, but we will provide directions for doing everything without one
that you can follow if you prefer.

## Getting the data

If you haven't already done so, clone the repository for this workshop:

```
git clone https://github.com/glaParaBio/helminth-apollo-workshop-2025
cd helminth-apollo-workshop-2025
```

You'll now have two data directories called `data/` and `jbrowse_data/` and the
docker configuration file `compose.yml`.

## Running Apollo

Now in the terminal, run `docker compose up`. You should see a stream of logs
from the Docker containers. If you do, Apollo is now running! You can use
<kbd>Ctrl</kbd> + <kbd>C</kbd> in the terminal to stop Apollo when you are done.

Open your web browser and navigate to [http://localhost/](http://localhost/).
You should see the JBrowse front page and clicking on *Linear genome view* will
take you to a new session, now empty. Next step is to [load some data](02-loading-data.md).

---

> [!NOTE]
> In case you need to remove the docker components (and data!) associated to Apollo:

```sh
cont_id=`docker ps -a | grep -P -w 'apollo-local-testing-client|ghcr.io/gmod/apollo-collaboration-server|mongo:7' | cut -f 1 -d ' '`
docker stop $cont_id
docker system prune
docker volume rm $(docker volume ls -qf dangling=true)
image_id=`docker images -a | grep -P -w 'apollo-local-testing-client|mongo|ghcr.io/gmod/apollo-collaboration-server' | awk '{print $3}'`
docker rmi -f $image_id
```
