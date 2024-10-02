/*
        scale-free network generator adapted from:
        https://github.com/visjs/vis-network/blob/master/examples/network/exampleUtil.js
        Copyright (c) 2014-2017 Almende B.V.
*/
function getScaleFreeNetwork(nodeCount) {
    const nodes = [];
    const edges = [];
    const connectionCount = [];
    // randomly create some nodes and edges
    for (let i = 0; i < nodeCount; i++) {
        nodes.push({
            id: i,
            label: String(i),
        });

        connectionCount[i] = 0;
        // create edges in a scale-free-network way
        if (i === 1) {
            const from = i;
            const to = 0;
            edges.push({
              from: from,
              to: to,
            });
            connectionCount[from]++;
            connectionCount[to]++;
        } else if (i > 1) {
            const conn = edges.length * 2;
            const rand = Math.floor(seededRandom() * conn);
            let cum = 0;
            let j = 0;
            while (j < connectionCount.length && cum < rand) {
                cum += connectionCount[j];
                j++;
            }

            const from = i;
            const to = j;
            edges.push({
                from: from,
                to: to,
            });
            connectionCount[from]++;
            connectionCount[to]++;
        }
    }
    return { nodes: nodes, edges: edges };
}

function sleep(ms) {
  return new Promise(resolve => setTimeout(resolve, ms));
}

async function animateMove(
    startX, startY, endX, endY,
    steps, interval, id
) {
    let xStep = (endX - startX) / steps;
    let yStep = (endY - startY) / steps;
    for (let i = 1; i <= steps; i++) {
        network.moveNode(
            id, startX + (i * xStep),
            startY + (i * yStep)
        );
        await sleep(interval);
    }
}

function randomMove() {
    let xRange = document.getElementById("network").clientWidth;
    let yRange = document.getElementById("network").clientHeight;
    let id = Math.floor(Math.random() * (nNodes + 1));
    let x, y;
    let xSign = Math.random() < 0.5 ? -1 : 1;
    let ySign = Math.random() < 0.5 ? -1 : 1
    x = Math.floor(Math.random() * xRange + xRange/2) * xSign;
    y = Math.floor(Math.random() * yRange + yRange/2) * ySign;
    let nodePosition = network.getPositions()[id];
    if (nodePosition !== null && nodePosition !== undefined) {
        animateMove(
            nodePosition.x, nodePosition.y, x, y,
            2000, 10, id
        );
    }
}

function scaleNetwork() {
    let networkHeight = screen.height - document.getElementById("welcome_text").getBoundingClientRect()['bottom'] - 1;
    document.getElementById("network").style.height = networkHeight + "px;"
    return networkHeight;
}