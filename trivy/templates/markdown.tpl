{{- if . }}
{{- range . }}
<h2>Trivy scan</h2>
{{- if (eq (len .Vulnerabilities) 0) }}
<h3>No Vulnerabilities found</h3>
{{- else }}
<h3>Vulnerabilities ({{ len .Vulnerabilities }})</h3>
<table>
    <tr>
        <th>Package</th>
        <th>ID</th>
        <th>Severity</th>
        <th>Installed Version</th>
        <th>Fixed Version</th>
    </tr>
    {{- range .Vulnerabilities }}
    <tr>
        <td><code>{{ escapeXML .PkgName }}</code></td>
        <td>{{ escapeXML .VulnerabilityID }}</td>
        <td>{{ escapeXML .Severity }}</td>
        <td>{{ escapeXML .InstalledVersion }}</td>
        <td>{{ escapeXML .FixedVersion }}</td>
    </tr>
    {{- end }}
</table>
{{- end }}
{{- if (eq (len .Misconfigurations ) 0) }}
<h3>No Misconfigurations found</h3>
{{- else }}
<h3>Misconfigurations</h3>
<table>
    <tr>
        <th>Type</th>
        <th>ID</th>
        <th>Check</th>
        <th>Severity</th>
        <th>Message</th>
    </tr>
    {{- range .Misconfigurations }}
    <tr>
        <td>{{ escapeXML .Type }}</td>
        <td>{{ escapeXML .ID }}</td>
        <td>{{ escapeXML .Title }}</td>
        <td>{{ escapeXML .Severity }}</td>
        <td>
          {{ escapeXML .Message }}
          <br><a href={{ escapeXML .PrimaryURL | printf "%q" }}>{{ escapeXML .PrimaryURL }}</a></br>
        </td>
    </tr>
    {{- end }}
</table>
{{- end }}
{{- end }}
{{- else }}
<h2>Trivy Returned Empty Report</h2>
{{- end }}
