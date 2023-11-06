#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include<stdbool.h>
#include <SDL.h>
#define WINDOW_WIDTH 1870      // largeur de la fenêtre pour afficher le piano
#define WINDOW_HEIGHT 400		// hauteur de la fenêtre pour afficher le piano
	//gcc main.c -o prog $(sdl2-config --cflags --libs)


#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])
#define K 1.0234


int32_t* wav2table(char* nom, int32_t* T_ampl,long int* pfreq, unsigned int* pOctpSec, unsigned int* pBpSample, unsigned long int* ptaille, unsigned int* pn_Can){
    FILE* f = fopen(nom,"rb");
    unsigned char* header = malloc(44*sizeof(char));
    fread(header, sizeof(char), 44, f);
    *ptaille = header[4] + header[5] * 16*16 + header[6] * 16*16*16*16 + header[7] * 16*16*16*16*16*16 + 8; //On lit les données en petit-boutiste
    for (unsigned int i=0; i<44; i++){
    }
    *pOctpSec = header[28] + header[29] * 16*16 + header[30] * 16*16*16*16 + header[31] * 16*16*16*16*16*16;
    *pfreq = header[24] + header[25] * 16*16 + header[26] * 16*16*16*16 + header[27] * 16*16*16*16*16*16;
    *pBpSample = header[34] + header[35] * 16*16;
    *pn_Can = header[22] + header[23]*16*16;
    free(header);
    fclose(f);
    FILE* f2 = fopen(nom,"rb");
    header = malloc(4*sizeof(char));
    fread(header, sizeof(char), 4, f2);
    unsigned char a[1] = {header[3]};
    unsigned char b = header[2];
    unsigned char c = header[1];
    unsigned char d = header[0];
    long unsigned int i = 4;
    //On se repère grâce aux caractères "DATA"
    while ((i<*ptaille) && ((a[0] != 6*16+1) || (b != 7*16+4) || (c != 6*16+1) || (d != 6*16+4))){
        d = c;
        c = b;
        b = a[0];
        fread(a, sizeof(char), 1, f2);
        i++;
    }
    fread(header, sizeof(char), 4, f2);
    long unsigned int t_data = header[0] + header[1] *16*16 + header[2] *16*16*16*16 + header[3] *16*16*16*16*16*16;
    free(header);
    unsigned char* T_data = malloc(t_data*sizeof(unsigned char));
    unsigned long int m = fread(T_data, sizeof(char), t_data, f2);
    fclose(f2);
    T_ampl = malloc(t_data*sizeof(int32_t));
    if (*pBpSample == 16){ //On convertit les octets  dans le bon format
        for (unsigned long int i=0; i<t_data; i+=2){
            if (T_data[i+1]>=128){
                T_ampl[i] = -32768 + ((T_data[i+1] - 128) * 256 + T_data[i]);
            }
            else{
                T_ampl[i] = 256 * T_data[i+1] + T_data[i];
            }
        }
    }
    else if (*pBpSample == 8){
        for (unsigned long int i=0; i<t_data; i+=1){
            T_ampl[i] = T_data[i] -128;
        }
    }
    else if (*pBpSample == 24){
        for (unsigned long int i=0; i<t_data; i+=3){
            if (T_data[i+2]>=128){
                T_ampl[i] = -8388608 + ((T_data[i+2] - 128) * 65536 + T_data[i+1] * 256 + T_data[i]);
            }
            else{
                T_ampl[i] = 65536 * T_data[i+2] + 256 * T_data[i+1] + T_data[i];
            }
        }
    }
    else if (*pBpSample == 32){
        for (unsigned long int i=0; i<t_data; i+=4){
            if (T_data[i+3]>=128){
                T_ampl[i] = -2147483648 + ((T_data[i+3] - 128) * 16777216 + T_data[i+2] * 65536 + T_data[i+1] * 256 + T_data[i]);
            }
            else{
                T_ampl[i] = 16777216 * T_data[i+3] + 65536 * T_data[i+2] + 256 * T_data[i+1] + T_data[i];
            }
        }
    }
    *ptaille = t_data;
    return T_ampl;
}

int32_t hamming(int32_t a, long unsigned int i, long unsigned int taille){
	return a*(0.54 + 0.46*sin(M_PI*i*1.0/taille));
}

double* FFT(double* T, long unsigned int taille, double* pftaille){
	*pftaille = pow(2,floor(log2(taille)+1));
	double* Tc = malloc(sizeof(double)*(*pftaille)*2);
	for (long unsigned int i = 0; i < taille; i++)
	{
	REAL(Tc,i) = T[i];
	IMAG(Tc,i) = 0.0;
        }
       for (long unsigned int i = taille; i < *pftaille; i++)
         {
         REAL(Tc,i) = 0.0;
         IMAG(Tc,i) = 0.0;
         }
       gsl_fft_complex_radix2_forward(Tc, 1, *pftaille);
       return Tc;
}
long unsigned int prod_spec(double* T, unsigned int H, long unsigned int taille, double ftaille){
	long unsigned int imax = 0;
	double m = log(REAL(T,taille/2-1)*REAL(T,taille/2-1)+IMAG(T,taille/2-1)*IMAG(T,taille/2-1));
	for (long unsigned int i=0; i<ftaille/2; i+=1){ //On peut s'arrêter à la moitié du tableau car il est symétrique
		double p=1.0;
		for (unsigned int l=1; l<=H; l+=1){ //Calcul du maximum
			if (i*l<ftaille/2){
				p*=REAL(T,i*l)*REAL(T,i*l)+IMAG(T,i*l)*IMAG(T,i*l);
			}
		}
		if (log(p)>m){
			m = log(p);
			imax = i;
		}
	}
	return imax;
}
long unsigned int f_max(double* T, long unsigned int taille, double* pftaille){
	double* TF = FFT(T, taille, pftaille);
	double if_max = prod_spec(TF, 3, taille, *pftaille);//Indice de fréquence max du produit spectral
	free(TF);
	return if_max;
}

int* Pierre(int* nombre_de_notes, char* nom){
    int32_t* T_ampl; //Tableau des amplitudes en int_32
    long int freq; //Fréquence d'échantillonnage
    unsigned int OctpSec; //Octets par secondes
    unsigned int BpSample; //Bits par échantillon
    unsigned long int taille; //Taille de la partie audio du fichier en octets
    unsigned int n_Can; //nombre de canaux
    unsigned long int* ptaille = &taille; //pointeur de la taille
    double ftaille = 0.0;
    double* pftaille = &ftaille; //Pointeur de la taille des tableaux de fréquences (on remplit les tableaux de 0 pour avoir une puissance de 2 pour la transformée de Fourier
    double T_notes_freq[88] = {27.5,29.1,30.9,32.7,34.6,36.7,28.9,41.2,43.7,46.2,29,51.9,55,58.3,61.7,65.4,69.3,73.4,77.8,82.4,87.3,92.5,98,103.8,110,116.5,123.5,130.8,138.6,146.8,155.6,164.8,174.6,185,196,207.7,220,233.1,247,261.6,277.2,293.7,311.2,329.6,349.2,370,392,415.3,440,466.2,493.9,523.3,554.4,587.3,622.3,659.3,698.5,740,784,830.6,880,932.3,987.8,1047,1109,1175,1245,1319,1397,1480,1568,1661,1760,1865,1975,2093,2217,2349,2489,2637,2794,2960,3136,3322,3520,3729,3951,4186}; //fréquences des notes
    T_ampl = wav2table(nom, T_ampl, &freq, &OctpSec, &BpSample, &taille, &n_Can); //Tableau des amplitudes (entiers)
    double duree = (1.0*taille)/OctpSec; //durée du morceau
    taille = taille/BpSample*8/n_Can; //taille en nombre d'échantillons
    long unsigned int n_morceaux = duree/0.1; //nombre de morceaux depetite durée à analyser
    long unsigned int long_morceaux = taille/n_morceaux; //nombre d'échantillons dans les petits morceaux à analyser
    double* T = malloc(sizeof(double)*taille); //tableau des amplitudes comprises entre -1 et 1
    double maxi = 0;
    for (unsigned long int k=0; k<n_morceaux; k++){ //on calcule les amplitudes après avoir passé le fenêtrage
    	for (unsigned long int i=0; i<long_morceaux; i++){
    		T_ampl[k*long_morceaux+i] = hamming(T_ampl[k*long_morceaux+i], i, long_morceaux);
    	}
    }
    for (unsigned long int i=0; i<taille; i++){ //on calcule l'amplitude maxi
    	if (abs(T_ampl[i*n_Can*BpSample/8]) > maxi){
    		maxi = abs(T_ampl[i*n_Can*BpSample/8]);
    	}
    }
    for (unsigned long int i=0; i<taille; i++){ //on normalise le morceau
    	T[i] = T_ampl[i*n_Can*BpSample/8]/maxi;
    }
    int* Tableau_final = malloc(sizeof(int)*n_morceaux); //On initialise le tableau final compenrtant toutes les notes
    *nombre_de_notes = n_morceaux;
    for (unsigned long int i=0; i<n_morceaux; i++){ //on traite chaque petit morceau un par un
    	long unsigned int ifreq_max = f_max(T+i*long_morceaux, long_morceaux, pftaille); //Fréquence du maximum du produit spectral sur cette partie
    	int k = 0;
    	while ((k<88) && ((T_notes_freq[k]*K<ifreq_max*freq/(*pftaille)) || (T_notes_freq[k]>K*ifreq_max*freq/(*pftaille)))){ //Si la fréquence correspond à une note de la gamme
    		k++;
    	}
    	if (k<88){
    		Tableau_final[i] = k+1; //la touche jouée est la k+1ème (convention en commençant à 1 pour les touches de piano

    	}
    	else{
    		Tableau_final[i] = 89; //Aucune note jouée


    	}
    }
    return Tableau_final; //On renvoie le taleau contenant toutes les notes
}


















void SDL_ExitWithError(const char *message);
int touche_blanche(int touche);
void colorier_piano(SDL_Renderer *renderer, int touche,SDL_Texture *texture);


int main(int argc, char **argv){
    	char* nom = "gamme.wav";// Musique à lire
	SDL_Window *window=NULL;     // Création et initialisation de la fenêtre
	SDL_Renderer *renderer=NULL ;	// Création et initialisation du rendu
	int N;	// Nombre de notes du morceau musique
	int* Tableau_final = Pierre(&N, nom);
	//SDL_Renderer *renderer_tmp=NULL ;

	if(SDL_Init(SDL_INIT_VIDEO)!=0)     // On vérifie que l'on a bien pu initialiser la vidéo

	{
		SDL_ExitWithError("Initialisation SDL");
	}
	//Création de la fenêtre ainsi que du rendu
	if(SDL_CreateWindowAndRenderer(WINDOW_WIDTH,WINDOW_HEIGHT,0,&window,&renderer)!=0){
		SDL_ExitWithError("Création de la fenêtre et du rendu échouée");
	}
	//SDL_CreateRenderer(&window,-1,0);

	SDL_Surface *image=NULL;  // création et initialisation de la surface qui va accueillir l'image
	SDL_Texture *texture=NULL;		//création et initialisation de la texture

	image=SDL_LoadBMP("pianof.bmp");      // Chargement du piano sur la surface

	// On vérifie que l'on a pu charger l'image et le cas contraire on détruit la fenêtre et le rendu.
	if (image==NULL){
		SDL_DestroyRenderer(renderer);
		SDL_DestroyWindow(window);
		SDL_ExitWithError("Impossible de charger l'image");
	}

	texture=SDL_CreateTextureFromSurface(renderer,image); // On crée la texture à partir de l'image du tableau
	SDL_Texture *texture_tmp=SDL_CreateTextureFromSurface(renderer,image);


	//SDL_FreeSurface(image); // On libère la surface car elle ne nous est plus utile.

	// On vérifie qu'on a bien pu charger la texture et le cas contraire on ferme le rendu et la fenêtre.
	if (texture==NULL){
		SDL_DestroyRenderer(renderer);
		SDL_DestroyWindow(window);
		SDL_ExitWithError("Impossible de charger la texture");

	}
	SDL_Rect rectangle; //  Création d'un rectangle
	if(SDL_QueryTexture(texture, NULL, NULL ,&rectangle.w, &rectangle.h)!=0 ){
		SDL_DestroyRenderer(renderer);
		SDL_DestroyWindow(window);
		SDL_ExitWithError("Impossible de charger la texture");

	}
	rectangle.x=(WINDOW_WIDTH-rectangle.w)/2;    //On centre l'image du piano
	rectangle.y=(WINDOW_HEIGHT-rectangle.h)/2;


	if(SDL_RenderCopy(renderer,texture, NULL,&rectangle)!=0){
		SDL_DestroyRenderer(renderer);
		SDL_DestroyWindow(window);
		SDL_ExitWithError("Impossible de charger la texture");

	}
	for(int j=0;j<=N-1;j++){
		SDL_Texture *texture_tmp=SDL_CreateTextureFromSurface(renderer,image);

		int touche_tmp=Tableau_final[j];   // on récupère la valeur de la touche à afficher
		colorier_piano(renderer,touche_tmp,texture);
		SDL_RenderPresent(renderer);
		SDL_Delay(92);  // Temps  pour calculer entre chaque note
		rectangle.x=(WINDOW_WIDTH-rectangle.w)/2;    //On centre l'image du piano
		rectangle.y=(WINDOW_HEIGHT-rectangle.h)/2;
		SDL_QueryTexture(texture_tmp, NULL, NULL ,&rectangle.w, &rectangle.h);   // on crée une texture temporaire pour conserver le piano vierge d'une part et pouvoir dessiner sur un autre piano en parallèle
		SDL_RenderCopy(renderer,texture_tmp, NULL,&rectangle);




		;


	}

	SDL_FreeSurface(image); //On libère la surfaace
	SDL_RenderPresent(renderer);
	SDL_Delay(5000);     // Délai de 5s

	SDL_DestroyTexture(texture);// On détruit les structures
	SDL_DestroyTexture(texture_tmp);


	SDL_DestroyRenderer(renderer);  // On détruit le rendu
	SDL_DestroyWindow(window);      // On détruit la fenêtre
	SDL_Quit();			// On quitte la  SDL



	return EXIT_SUCCESS;
}

// Fonction permettant d'afficher simplement un message d'erreur et de quitter la SDL en cas de problèmes

void SDL_ExitWithError(const char *message){
	SDL_Log("ERREUR : %s > %s\n", message , SDL_GetError());
	SDL_Quit();
	exit(EXIT_FAILURE);

}


// Fonction renvoyant un booléen indiquant si la touche est noire ou non .
bool touche_noire(int touche){
	if ((touche==2)||(touche==5)||(touche==7)||(touche==10)||(touche==12)||(touche==14)||(touche==17)||(touche==19)||(touche==22)||(touche==24)||(touche==26)||(touche==29)||(touche==31)||(touche==34)||(touche==36)||(touche==38)||(touche==41)||(touche==43)||(touche==46)||(touche==48)||(touche==50)||(touche==53)||(touche==55)||(touche==58)||(touche==60)||(touche==62)||(touche==65)||(touche==67)||(touche==70)||(touche==72)||(touche==74)||(touche==77)||(touche==79)||(touche==82)||(touche==84)||(touche==86)){
		return true;
	}
	return false;
}

// Fonction qui permet de colorier en rouge les touches jouées sur le piano
void colorier_piano(SDL_Renderer *renderer, int touche,SDL_Texture *texture){
	//On regarde si la touche jouée est une touche noire
	SDL_SetRenderTarget(renderer,texture); //Permet de dessiner sur la texture et non le rendu
	// Lorsque la touche vaut 89 cela signifie qu'on ne joue pas de notes dans la musique, on colorie seulement un petit point en rouge
	if(touche==89){
	SDL_Color rouge={255,1,1,255};
	SDL_SetRenderDrawColor(renderer, rouge.r, rouge.g, rouge.b, rouge.a);
	SDL_RenderDrawPoint(renderer,0,0);}

	else if(touche_noire(touche)){
		SDL_Color rouge={255,1,1,255}; // Structure représentant la couleur rouge
		SDL_SetRenderDrawColor(renderer, rouge.r, rouge.g, rouge.b, rouge.a);     // Initialise la couleur à rouge pour le dessin sur le rendu
		SDL_Rect rect={21*touche+1,22,26,256};    // Défini le rectangle à colorier en rouge
		SDL_RenderFillRect(renderer,&rect); // Colorie en rouge le rectangle en question
	}
	//On fait la même chose pour une touche blanche
	else{
		int touche_b=touche_blanche(touche);
		SDL_Color rouge={255,1,1,255};
		SDL_SetRenderDrawColor(renderer, rouge.r, rouge.g, rouge.b, rouge.a);
		SDL_Rect rect={21.9+(35.1)*(touche_b-1),279,33.9,102};

		SDL_RenderFillRect(renderer,&rect);

	}
	SDL_SetRenderTarget(renderer,NULL);


}
// Fonction permettant de définir une valeur touche blanche comrprise entre 1 et 52 afin de dessiner plus facilement sur la texture
int touche_blanche(int touche){
	int touche_b;
	int touche_tmp=touche-3;      // Pour trouver l'octave dans laquelle on se trouve, on enlève les 3 premières touches afin d'avoir une octave pleine
	int octave=touche_tmp/12;     // Donne l'octave dans laquelle on se trouve
	int place=touche_tmp%12;      // Donne la place de la touche blanche dans son octave


	if (touche==1){
		touche_b=1;

	}

	else if(touche==3){
		touche_b=touche-1;


	}

	else if(place==1){
		touche_b=touche-1-5*octave;

	}

	else if(place==3){
		touche_b=touche-2-5*octave;

	}

	else if((place==5)||(place==6)){
		touche_b=touche-3-5*octave;

	}

	else if(place==8){
		touche_b=touche-4-5*octave;

	}

	else if(place==10){
		touche_b=touche-5-5*octave;

	}

	else if(place==0){
		touche_b=touche-5*(octave+1);

	}
	}







