<%@ page contentType="text/html;charset=big5"%>
<%@ page import="java.sql.*"%>
<%@ page import="java.io.*"%>
<%@ page import="java.util.*"%>
<%@ page import="java.text.*"%>

<head>

	<link rel="icon" href="http://syslab.nchu.edu.tw/IIIDB/images/syslab.ico" type="image/x-icon" />
    <link rel="shortcut icon" href="http://syslab.nchu.edu.tw/IIIDB/images/syslab.ico" type="image/x-icon" />
	<link rel="bookmark" href="http://syslab.nchu.edu.tw/IIIDB/images/syslab.ico" type="image/x-icon" />

    <title>IIIDB - welcome</title>
    
    <meta http-equiv="content-type" content="application/xhtml+xml; charset=UTF-8" />
    
        <link rel="stylesheet" type="text/css" href="css/simple.css" media="screen, tv, projection" title="Default" />
	<script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.4/jquery.min.js"></script>
	<script>
		!window.jQuery && document.write('<script src="/jquery/jquery-1.4.3.min.js"><\/script>');
	</script>
	<script type="text/javascript" src="./jquery/fancybox/jquery.mousewheel-3.0.4.pack.js"></script>
	<script type="text/javascript" src="./jquery/fancybox/jquery.fancybox-1.3.4.pack.js"></script>
	<link rel="stylesheet" type="text/css" href="./jquery/fancybox/jquery.fancybox-1.3.4.css" media="screen" />

	<script type="text/javascript">
		$(document).ready(function() {
			$("a[rel=example_group]").fancybox({
				'transitionIn'		: 'none',
				'transitionOut'		: 'none',
				'titlePosition' 	: 'over',
				'titleFormat'		: function(title, currentArray, currentIndex, currentOpts) {
					return '<span id="fancybox-title-over">Image ' + (currentIndex + 1) + ' / ' + currentArray.length + (title.length ? ' &nbsp; ' + title : '') + '</span>';
				}
			});		
		});
	</script>
    

<script type="text/javascript" src="jquery/jquery.tablesorter.js"></script>
<script type="text/javascript">
	$(function() {
		$("table").tablesorter({widgets: ['zebra']});
	});
</script>
<style type="text/css">
body {
	background-image: url(./images/bg1.jpg);
	background-repeat: no-repeat;
}
</style>
</head>
<body id="your-site-id">

<div id="page" class="classname">
    
    <div id="header">
        <h1><a href="index.html">IIIDB</a></h1>
        <strong>A database for Isoform-Isoform Interactions and isoform network modules.</strong></div>
    
    <div id="wrapper"> 
        <div id="right">
        <div id="content">
        
            <div id="path">
                You are here: <a href="index.html">Home</a> 
                &raquo; <a href="Enrich.html">Module Search</a>
            </div>
    
            <div id="main">

                <h1>Module Search</h1>
                <br />




<%
request.setCharacterEncoding("big5");
String ID = request.getParameter("TranGeneID");
ID = ID.trim();
Class.forName("com.mysql.jdbc.Driver");
Connection conn = DriverManager.getConnection("jdbc:mysql://140.120.203.50:3306/ting?useUnicode=true&characterEncoding=big5","ting", "ting2010");
Statement stmt = conn.createStatement();

NumberFormat f = NumberFormat.getNumberInstance();
f.setMaximumFractionDigits(5);

Set<String> modid_set = new HashSet<String>();
String genesym = "null";
String geneid = "null";
String get_mod_sql="select modid, genesym, geneid from isoform_mod WHERE geneid= '" + ID + "' OR genesym= '" + ID + "';" ;
if (stmt.execute(get_mod_sql))   {
   ResultSet rs = stmt.getResultSet();
   while (rs.next())
    {
		modid_set.add(rs.getString("modid"));
		genesym = rs.getString("genesym");
		geneid = rs.getString("geneid");
	}
   rs.close();
}

if (genesym.equals("null") && geneid.equals("null"))
{
%>	
	<h3>There are 0 module in your request: <%=ID%></h3>
<%
}else{
%>
	<h3>There are <%out.print(modid_set.size());%> modules in your request: 
	<a href="http://www.ncbi.nlm.nih.gov/gene?term=<%=ID%>" target=_blank><%=genesym%> (GeneID: <%=geneid%>)</a>
	</h3>
<%
}
for (String m: modid_set){
%>
<br />
<h3>Module ID: <%=m%></h3>
<a rel="example_group" href="./images/module/<%=m%>_B.gif" title="<%=m%>"><img alt="" src="./images/module/<%=m%>_S.gif" /></a> 

<%
	NumberFormat nf1 = NumberFormat.getInstance();
    nf1.setMaximumFractionDigits(2);
   	nf1.setMinimumFractionDigits(2);
	File ff = new File("../webapps/IIIDB/images/module/", m+"_B.pdf");
	double size = 0;
	if(ff.exists()){
		size = ff.length()/1024.0;
	}
%>	
<a href="./images/module/<%=m%>_B.pdf" title="<%=m%>" target=_blank><img src="images/zoom_in.png" width="16" height="16">Full figure (<%= nf1.format(size)%>KB)</a>

<%
String sql2="select * from mod_complex_annot WHERE modid= '" + m + "' order by pValue;" ;
if (stmt.execute(sql2))   {
	ResultSet rs = stmt.getResultSet();
	if (rs.isAfterLast()!=rs.isBeforeFirst()){
		
%>


<br />
<h3>Complex Enrichment Analysis: </h3>

<table cellspacing="0" class="tablesorter">
<thead>
    <TR>
    <th width="12%" nowrap="nowrap"><B>Module ID</B></B></th>
    <th><b>Complex</b></th>
    <th><B>Associated genes</B></th>
    <th width="10%" nowrap="nowrap"><B>P-value</B></B></th>
    </TR>   
</thead>
<tbody>
    <%
	NumberFormat nf = NumberFormat.getInstance();
    nf.setMaximumFractionDigits(8);//若小數點超過四位，則第五位~四捨五入
   	nf.setMinimumFractionDigits(8);//若小數點不足二位，則補足二位
    while (rs.next())
    {
		%>
        <TR>
        <TD nowrap="nowrap"><B><%= rs.getString("modid") %></B></TD>
        <TD><B><%= rs.getString("complex") %></B></TD>
        <TD><B><%
		String associGene = rs.getString("associGene");
		String[] items = associGene.split(", ");
		for (int i = 0; i < items.length; i++){
			%>
			<a href="Enrich.jsp?TranGeneID=<%= items[i] %>"><%= items[i] %></a>
            <%
		out.print(" ");
        }
		%></B></TD>
		<%
		double pvalue = rs.getDouble("pValue");
		if (pvalue >= 0.00000001){
			//pvalue = nf.format(pvalue);
			%>
        <TD width="10%" nowrap="nowrap"><B><%= nf.format(pvalue) %></B></TD> 
       		<%
		}else{
			//pvalue = nf.format(pvalue);
			%>
        <TD width="10%" nowrap="nowrap"><B> &lt;1e-8 </B></TD> 
       		<%
		}
		%>
		</TR>
        <%
    }
   
	rs.close();
	
%>
</tbody>
    </TABLE>
    
<%
}
}
%>  

<p>

<%
String sql="select * from mod_path_annot WHERE modid= '" + m + "' order by pValue;" ;
if (stmt.execute(sql))   {
	ResultSet rs = stmt.getResultSet();
	if (rs.isAfterLast()!=rs.isBeforeFirst()){
		
%>

<h3>Pathway Enrichment Analysis: </h3>

<table cellspacing="0" class="tablesorter">
<thead>
    <TR>
    <th width="12%" nowrap="nowrap"><B>Module ID</B></B></th>
    <th><b>Pathway</b></th>
    <th><B>Associated genes</B></th>
    <th width="10%" nowrap="nowrap"><B>P-value</B></B></th>
    </TR>   
</thead>
<tbody>
    <%
	NumberFormat nf = NumberFormat.getInstance();
    nf.setMaximumFractionDigits(8);//若小數點超過四位，則第五位~四捨五入
   	nf.setMinimumFractionDigits(8);//若小數點不足二位，則補足二位
    while (rs.next())
    {
		%>
        <TR>
        <TD nowrap="nowrap"><B><%= rs.getString("modid") %></B></TD>
        <TD><B><%= rs.getString("pathway") %></B></TD>
        <TD><B><%
		String associGene = rs.getString("associGene");
		String[] items = associGene.split(", ");
		for (int i = 0; i < items.length; i++){
			%>
			<a href="Enrich.jsp?TranGeneID=<%= items[i] %>"><%= items[i] %></a>
            <%
		out.print(" ");
        }
		%></B></TD>
		<%
		double pvalue = rs.getDouble("pValue");
		if (pvalue >= 0.00000001){
			//pvalue = nf.format(pvalue);
			%>
        <TD width="10%" nowrap="nowrap"><B><%= nf.format(pvalue) %></B></TD> 
       		<%
		}else{
			//pvalue = nf.format(pvalue);
			%>
        <TD width="10%" nowrap="nowrap"><B> &lt;1e-8 </B></TD> 
       		<%
		}
		%>
		</TR>
        <%
    }
   
	rs.close();
	
%>
</tbody>
    </TABLE>
    
<%
}
}
%>  

<p>

  <%
String sql1="select * from mod_GO_annot WHERE modid= '" + m + "' order by pValue;" ;

if (stmt.execute(sql1))   {
   ResultSet rs1 = stmt.getResultSet();
	if (rs1.isAfterLast()!=rs1.isBeforeFirst()){
%>

<h3>GO Enrichment Analysis: </h3>

<table  cellspacing="0" class="tablesorter">
    <thead>
    <TR>
    <th width="12%" nowrap="nowrap"><B>Module ID</B></th>
    <th><B>GO</B></th>
    <th><B>Associated genes</B></th>
    <th width="10%" nowrap="nowrap"><B>P-value</B></B></th>
    </TR>   
    </thead>
    <tbody>
    <%
    NumberFormat nf = NumberFormat.getInstance();
    nf.setMaximumFractionDigits(8);//若小數點超過四位，則第五位~四捨五入
   	nf.setMinimumFractionDigits(8);//若小數點不足二位，則補足二位
    while (rs1.next())
    {
        %>
        <TR>
        <TD nowrap="nowrap"><B><%= rs1.getString("modid") %></B></TD>
        <TD><B><%= rs1.getString("GO") %></B></TD>
        <TD><B><%
		String associGene = rs1.getString("associGene");
		String[] items = associGene.split(", ");
		for (int i = 0; i < items.length; i++){
			%>
			<a href="Enrich.jsp?TranGeneID=<%= items[i] %>"><%= items[i] %></a>
            <%
		out.print(" ");
        }
		%></B></TD>
        <%
		double pvalue = rs1.getDouble("pValue");
		if (pvalue >= 0.00000001){
			//pvalue = nf.format(pvalue);
			%>
        <TD width="10%" nowrap="nowrap"><B><%= nf.format(pvalue) %></B></TD> 
       		<%
		}else{
			//pvalue = nf.format(pvalue);
			%>
        <TD width="10%" nowrap="nowrap"><B>&lt;1e-8 </B></TD> 
       		<%
		}
		%>
		</TR>
        <%
    }
}
	rs1.close();      
%>
</tbody>
    </TABLE>
    
<%


}
}
    
    stmt.close();  
    conn.close();  
    %>  

<br>

<br>
                
               
                
            </div>

        </div>
        </div>
        <div id="left">
          <div id="nav">
            <h3>Site Map</h3>
            <div class="inner_box">
              <ul id="mb73g0ebul_table" class="mb73g0ebul_menulist" style="width: 180px; height: 150px;">
                <li class="spaced_li"><a href="index.html#index"><img id="mbi_mb73g0_1" src="menu_images/mb_home.gif" name="mbi_mb73g0_1" width="160" height="31" style="vertical-align: bottom;" border="0" alt="Home" title="" /></a></li>
                <li class="spaced_li"><a href="Interact.html"><img id="mbi_mb73g0_2" src="menu_images/mb_interaction_search.gif" name="mbi_mb73g0_2" width="160" height="31" style="vertical-align: bottom;" border="0" alt="Interaction Search" title="" /></a></li>
                <li class="spaced_li"><a href="Enrich.html"><img id="mbi_mb73g0_3" src="menu_images/mb_module_search.gif" name="mbi_mb73g0_3" width="160" height="31" style="vertical-align: bottom;" border="0" alt="Module Search" title="" /></a></li>
                <li class="spaced_li"><a href="index.html#about"><img id="mbi_mb73g0_5" src="menu_images/mb_about_us.gif" name="mbi_mb73g0_5" width="160" height="31" style="vertical-align: bottom;" border="0" alt="About Us" title="" /></a></li>
                <li class="spaced_li"><a href="index.html#tutorial"><img id="mbi_mb73g0_4" src="menu_images/mb_tutorial.gif" name="mbi_mb73g0_4" width="160" height="31" style="vertical-align: bottom;" border="0" alt="Tutorial" title="" /></a></li>
              </ul>
              <p>&nbsp;</p>
              <p>&nbsp;</p>
            </div>
            <link rel="stylesheet" href="menu_images/mbcsmb73g0.css" type="text/css" />
            <script type="text/javascript" src="menu_images/mbjsmb73g0.js"></script>
            <!-- i love clean source code -->
          </div>
          <div class="left_box">
            <h3>News and Updates</h3>
            <div class="inner_box">
              <p> <b>July 9, 2011:</b><br />
                New releases and related tools <br />
                <br />
                <a href="index.html">Read more...</a> </p>
            </div>
          </div>
          <div class="left_box">
            <h3>Related Tools</h3>
            <div class="inner_box">
              <p><a href="http://jigsaw.w3.org/css-validator/check/referer"> <img style="border:0;width:88px;height:31px"
        src="http://jigsaw.w3.org/css-validator/images/vcss-blue"
        alt="Valid CSS!" /> </a> </p>
<p> <a href="http://validator.w3.org/check?uri=referer"><img
      src="http://www.w3.org/Icons/valid-xhtml10" alt="Valid XHTML 1.0 Strict" height="31" width="88" /></a> </p>
            </div>
          </div>
        </div>
    </div>
    
    <div id="footer">
        <p>
        Copyright &copy; 2011 <a href="http://syslab.nchu.edu.tw/" target=_blank>Syslab</a>.</p>
    </div>

</div>
<a style="display:scroll;position:fixed;bottom:0px;right:5px;" href="#" title="" onFocus="if(this.blur)this.blur()"><img alt='' border='0' align="right" onMouseOver="this.src='http://4.bp.blogspot.com/-AKP9Gl_ets0/TmzrRBJPg_I/AAAAAAAADpY/M3KVvsDNJbA/s1600/Top.png'" src="http://1.bp.blogspot.com/-_gUa3K1wDNo/TmzrHROhBVI/AAAAAAAADpU/vOKm_7zL8DQ/s1600/Top_medium.png" onMouseOut="this.src='http://1.bp.blogspot.com/-_gUa3K1wDNo/TmzrHROhBVI/AAAAAAAADpU/vOKm_7zL8DQ/s1600/Top_medium.png'" /></a><!--BACKtoTOP-START--><!--BACKtoTOP-STOP-->
</body>
</html>