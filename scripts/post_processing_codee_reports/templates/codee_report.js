const chartData = ${CHART_DATA};

function getFileUrl(path) {
  if (!path) return null;
  const normalizedPath = path.replace(/\\/g, '/');
  return 'file://' + normalizedPath;
}

function initLinksTable() {
  const linksSection = document.getElementById('linksSection');
  const linksTableBody = document.getElementById('linksTableBody');

  if (chartData.hasLinks && chartData.reportLinks) {
    chartData.reportLinks.forEach((linkData, index) => {
      const date = chartData.labels[index];
      const row = document.createElement('tr');

      if (linkData && linkData.link) {
        const fileUrl = getFileUrl(linkData.link);
        const typeLabel = linkData.type === 'json' ? 'JSON' : 'HTML';
        const fileName = linkData.link.split('/').pop() || linkData.link;
        row.innerHTML = `
          <td>${date}</td>
          <td>Report ${index + 1} (${typeLabel})</td>
          <td><a href="#" onclick="window.open('${fileUrl}', '_blank'); return false;">${fileName}</a></td>
        `;
      } else {
        row.innerHTML = `
          <td>${date}</td>
          <td>Report ${index + 1}</td>
          <td><span class="no-links">No link available</span></td>
        `;
      }
      linksTableBody.appendChild(row);
    });
  } else {
    linksSection.style.display = 'none';
  }
}

function getLinkAtIndex(idx) {
  if (chartData.reportLinks && chartData.reportLinks[idx]) {
    return chartData.reportLinks[idx].link;
  }
  return null;
}

function getLinkTypeAtIndex(idx) {
  if (chartData.reportLinks && chartData.reportLinks[idx]) {
    return chartData.reportLinks[idx].type || 'unknown';
  }
  return null;
}

function createTotalChart() {
  const ctx = document.getElementById('totalChart').getContext('2d');
  new Chart(ctx, {
    type: 'line',
    data: {
      labels: chartData.labels,
      datasets: [{
        label: 'Total Findings',
        data: chartData.total,
        borderColor: '#667eea',
        backgroundColor: 'rgba(102, 126, 234, 0.1)',
        fill: true,
        tension: 0.3,
        pointRadius: 6,
        pointHoverRadius: 8,
        borderWidth: 3
      }]
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      plugins: {
        legend: { display: false },
        tooltip: {
          backgroundColor: 'rgba(0,0,0,0.8)',
          padding: 12,
          cornerRadius: 8,
          callbacks: {
            afterLabel: function(context) {
              const idx = context.dataIndex;
              const link = getLinkAtIndex(idx);
              if (link) {
                const type = getLinkTypeAtIndex(idx);
                return `Click to open ${type} report`;
              }
              return '';
            }
          }
        }
      },
      onClick: function(event, elements) {
        if (elements.length > 0) {
          const idx = elements[0].index;
          const link = getLinkAtIndex(idx);
          if (link) {
            const fileUrl = getFileUrl(link);
            window.open(fileUrl, '_blank');
          }
        }
      },
      scales: {
        y: {
          beginAtZero: true,
          grid: { color: 'rgba(0,0,0,0.05)' }
        },
        x: {
          grid: { display: false }
        }
      }
    }
  });
}

function createCheckerChart() {
  const checkerPalette = [
    '#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4',
    '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff',
    '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
    '#000075', '#a9a9a9', '#ffffff', '#000000', '#a93226', '#2471a3'
  ];

  const checkerDatasets = Object.entries(chartData.checkers).map(([checker, data], i) => {
    const color = checkerPalette[i % checkerPalette.length];
    return {
      label: checker,
      data: data,
      borderColor: color,
      backgroundColor: color,
      tension: 0.3,
      pointRadius: 4,
      pointHoverRadius: 6,
      fill: false
    };
  });

  const ctx = document.getElementById('checkerChart').getContext('2d');
  new Chart(ctx, {
    type: 'line',
    data: { labels: chartData.labels, datasets: checkerDatasets },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      interaction: {
        mode: 'index',
        intersect: false
      },
      plugins: {
        legend: {
          position: 'right',
          labels: { boxWidth: 12, padding: 15 }
        },
        tooltip: {
          backgroundColor: 'rgba(0,0,0,0.8)',
          padding: 12,
          cornerRadius: 8
        }
      },
      scales: {
        y: {
          beginAtZero: true,
          grid: { color: 'rgba(0,0,0,0.05)' }
        },
        x: {
          grid: { display: false }
        }
      }
    }
  });
}

function createPriorityChart() {
  const priorityDatasets = Object.entries(chartData.priorities)
    .filter(([_, data]) => data.some(v => v > 0))
    .map(([priority, data]) => ({
      label: priority,
      data: data,
      backgroundColor: chartData.priorityColors[priority] + 'CC',
      borderColor: chartData.priorityColors[priority],
      borderWidth: 1
    }));

  const ctx = document.getElementById('priorityChart').getContext('2d');
  new Chart(ctx, {
    type: 'line',
    data: { labels: chartData.labels, datasets: priorityDatasets },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      interaction: {
        mode: 'index',
        intersect: false
      },
      plugins: {
        legend: {
          position: 'top',
          labels: { boxWidth: 12, padding: 15 }
        },
        tooltip: {
          backgroundColor: 'rgba(0,0,0,0.8)',
          padding: 12,
          cornerRadius: 8
        }
      },
      scales: {
        y: {
          stacked: false,
          beginAtZero: true,
          grid: { color: 'rgba(0,0,0,0.05)' }
        },
        x: {
          grid: { display: false }
        }
      }
    }
  });
}

document.addEventListener('DOMContentLoaded', function() {
  initLinksTable();
  createTotalChart();
  createCheckerChart();
  createPriorityChart();
});
