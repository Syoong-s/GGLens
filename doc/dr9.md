Galaxy catalogs:

DESIDR9_galaxy: galaxy catalog 

        Number: 138315649 galaxies(~18000 deg^2) totally including 2.7% specz(~3700deg^2).
        no redshift cut.
        data format:--------------------------------------------------------------------------------------------------------------
            column 0: galaxy ID
                      NGC: 0-31450383 in north_ngc
                           31450384-72687991 in south_ngc
                      SGC: 72687992-138315648 in south_sgc
            column 1: ra
            column 2: dec
            column 3: redshift
            column 4: z-band absolute magnitude  (with k-correction_{0.5}^{z})
            column 5: z-band luminosity  (unit 10^10Lsun/h/h, with k-correction_{0.5}^{z}, Mr_{sun}=4.83)
            column 6: r-band absolute magnitude  (with k-correction_{0.5}^{r})
            column 7: r-band luminosity  (unit 10^10Lsun/h/h, with k-correction_{0.5}^{r}, Mr_{sun}=5.39)
            column 8: redshift error  (spec-z 0.0001; no redshift -99; photo-z other value.)
            column 9: log stellar mass
            column 10: z-band K-correction  (to z=0.5)
            column 11: z-band apparent magnitude
            column 12: r-band K-correction  (to z=0.5)
            column 13: r-band apparent magnitude
            column 14: select (-1:706288galaxies; 0:11020405galaxies; 1:122976436galaxies; 2:2213848galaxies; 3:1398672galaxies)
                        [0] Objects in NGC's two tails, summaried as 'select' >= -1.
                        [1] Objects with 5-band and reliable detection, summaried as 'select' >= 1. 
                            FLUX_[g,r,z,w1,w2] > 0,  FRACFLUX_[g,r,z]   < 0.5, 
                            FRAC_IN_[g, r, z]  > 0.3, FRACMASKED_[g,r,z] < 0.4.
                        [2] Objects with spec redshift [SV3+main+archives], summaried as 'select' >= 2
                        [3] Objects with DESI spec redshift [SV3+main], summaried as 'select' >= 3

----------------------------------------------

Group catalogs:

DESIDR9_NGC_group: group catalog at NGC
DESIDR9_SGC_group: group catalog at SGC

data format:

column 1: group ID
column 2: richness
column 3: ra
column 4: dec
column 5: redshift
column 6: log M_h/(M_sun/h)
column 7: L_group (10^10Lsun/h/h) 

----------------------------------------------

member galaxy files:

1. The galaxies are listed in `iDESIDR9_?GC_1' with format:
     1. galaxy id;  (galaxy ID in DESIDR9_galaxy)
     2. group id; 
     3. brightest galaxy = 1; other galaxies in groups = 2 - z-band
     4. auxiliary ID, not used.
   
