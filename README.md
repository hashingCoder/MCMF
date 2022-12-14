# MCMF

This is the sourse code of MCMF, as described in our paper:
 
[Learning to Hash Information Network Through MaxCut Formulation](https://)

## Data Sets


<table>
    <tr>
    	<th>Data set</th>
        <th>#Nodes</th>
        <th>#Edges</th>
        <th>#Classes</th>
    </tr>
    <tr>
    	<td> <a href="https://cloud.tsinghua.edu.cn/d/fd3717d9ee78440e800f/"> Protein-Protein Interaction </a> </td>
        <td> 3,890 </td>
        <td> 76,584 </td>
        <td> 50 </td>
    </tr>
    <tr>
    	<td> <a href="https://cloud.tsinghua.edu.cn/d/cb62b5b4224a4de08a02/"> BlogCatalog </a> </td>
        <td> 10,312 </td>
        <td> 333,983 </td>
        <td> 40 </td>
    </tr>
    <tr>
    	<td> <a href="https://cloud.tsinghua.edu.cn/d/a26619b0b45e4d1181c9/"> Wikipedia </a>  </td>
        <td> 4.777 </td>
        <td> 184,812 </td>
        <td> 39 </td>
    </tr>
    <tr>
    	<td> <a href="https://cloud.tsinghua.edu.cn/d/863da94f520844cbab90/"> Flickr </a> </td>
        <td> 80,513 </td>
        <td> 5,899,882 </td>
        <td> 195 </td>
    </tr>
</table>


## Usage
1. Download Data Sets 
2. Get Matrix Factorization S by [NetMF](https://github.com/xptree/NetMF) `./netmf.py` using function ` direct_compute_deepwalk_matrix` for T=1 and `approximate_deepwalk_matrix` for T=10, then save the obtain matrix S to `./data/network` 
3. Change the parameters in `./main_script.m`, i.e., dataset, gamma and dim
4. Run code in `./main_script.m`

## Cite

Please cite our paper if you use this code in your own work:

```

```