  "  ]   k820309    w          19.1        ���c                                                                                                          
       mod_tvd.f90 MOD_TVD          @       �                                  
                                                                                                             %         @                                                           #LHS    #RHS              
                                                     #MPI_INFO              
                                                     #MPI_INFO    %         @                                                           #LHS    #RHS 
             
                                                     #MPI_REQUEST 	             
                                  
                   #MPI_REQUEST 	   %         @                                                           #LHS    #RHS              
                                                     #MPI_COMM              
                                                     #MPI_COMM    %         @                                                           #LHS    #RHS              
                                                     #MPI_WIN              
                                                     #MPI_WIN    %         @                                                           #LHS    #RHS              
                                                     #MPI_FILE              
                                                     #MPI_FILE    %         @                                                           #LHS    #RHS              
                                                     #MPI_MESSAGE              
                                                     #MPI_MESSAGE    %         @                                                           #LHS    #RHS              
                                                     #MPI_ERRHANDLER              
                                                     #MPI_ERRHANDLER    %         @                                                           #LHS     #RHS "             
                                                      #MPI_GROUP !             
                                  "                   #MPI_GROUP !   %         @                                #                           #LHS $   #RHS &             
                                  $                   #MPI_DATATYPE %             
                                  &                   #MPI_DATATYPE %   %         @                                '                           #LHS (   #RHS *             
                                  (                   #MPI_OP )             
                                  *                   #MPI_OP )   %         @                                +                           #LHS ,   #RHS -             
                                  ,                   #MPI_INFO              
                                  -                   #MPI_INFO    %         @                                .                           #LHS /   #RHS 0             
                                  /                   #MPI_REQUEST 	             
                                  0                   #MPI_REQUEST 	   %         @                                1                           #LHS 2   #RHS 3             
                                  2                   #MPI_COMM              
                                  3                   #MPI_COMM    %         @                                4                           #LHS 5   #RHS 6             
                                  5                   #MPI_WIN              
                                  6                   #MPI_WIN    %         @                                7                           #LHS 8   #RHS 9             
                                  8                   #MPI_FILE              
                                  9                   #MPI_FILE    %         @                                :                           #LHS ;   #RHS <             
                                  ;                   #MPI_MESSAGE              
                                  <                   #MPI_MESSAGE    %         @                                =                           #LHS >   #RHS ?             
                                  >                   #MPI_ERRHANDLER              
                                  ?                   #MPI_ERRHANDLER    %         @                                @                           #LHS A   #RHS B             
                                  A                   #MPI_GROUP !             
                                  B                   #MPI_GROUP !   %         @                                C                           #LHS D   #RHS E             
                                  D                   #MPI_DATATYPE %             
                                  E                   #MPI_DATATYPE %   %         @                                F                           #LHS G   #RHS H             
                                  G                   #MPI_OP )             
                                  H                   #MPI_OP )                                            I                   
                &                                                    @                               J                   
                &                                                    @                               K                   
                &                                                    @                               L                   
                &                                                    @                               M                   
                &                                                    @                               N                   
                &                                                    @                               O                   
                &                                           #         @                                   P                    #SETUP_TVD%NCV Q   #SETUP_TVD%M R                                                                                    Q                                                       R                             @                           )     '                    #MPI_VAL S                �                               S                                    @                           %     '                    #MPI_VAL T                �                               T                                    @                           !     '                    #MPI_VAL U                �                               U                                    @                                '                    #MPI_VAL V                �                               V                                    @                                '                    #MPI_VAL W                �                               W                                    @                                '                    #MPI_VAL X                �                               X                                    @                                '                    #MPI_VAL Y                �                               Y                                    @                                '                    #MPI_VAL Z                �                               Z                                    @                           	     '                    #MPI_VAL [                �                               [                                    @                                '                    #MPI_VAL \                �                               \                      �         fn#fn    �   @   J   MOD_PREC    �   p       SP+MOD_PREC %   l  b       INFOEQ+MPI_CONSTANTS )   �  V   a   INFOEQ%LHS+MPI_CONSTANTS )   $  V   a   INFOEQ%RHS+MPI_CONSTANTS (   z  b       REQUESTEQ+MPI_CONSTANTS ,   �  Y   a   REQUESTEQ%LHS+MPI_CONSTANTS ,   5  Y   a   REQUESTEQ%RHS+MPI_CONSTANTS %   �  b       COMMEQ+MPI_CONSTANTS )   �  V   a   COMMEQ%LHS+MPI_CONSTANTS )   F  V   a   COMMEQ%RHS+MPI_CONSTANTS $   �  b       WINEQ+MPI_CONSTANTS (   �  U   a   WINEQ%LHS+MPI_CONSTANTS (   S  U   a   WINEQ%RHS+MPI_CONSTANTS %   �  b       FILEEQ+MPI_CONSTANTS )   
  V   a   FILEEQ%LHS+MPI_CONSTANTS )   `  V   a   FILEEQ%RHS+MPI_CONSTANTS (   �  b       MESSAGEEQ+MPI_CONSTANTS ,     Y   a   MESSAGEEQ%LHS+MPI_CONSTANTS ,   q  Y   a   MESSAGEEQ%RHS+MPI_CONSTANTS +   �  b       ERRHANDLEREQ+MPI_CONSTANTS /   ,  \   a   ERRHANDLEREQ%LHS+MPI_CONSTANTS /   �  \   a   ERRHANDLEREQ%RHS+MPI_CONSTANTS &   �  b       GROUPEQ+MPI_CONSTANTS *   F	  W   a   GROUPEQ%LHS+MPI_CONSTANTS *   �	  W   a   GROUPEQ%RHS+MPI_CONSTANTS )   �	  b       DATATYPEEQ+MPI_CONSTANTS -   V
  Z   a   DATATYPEEQ%LHS+MPI_CONSTANTS -   �
  Z   a   DATATYPEEQ%RHS+MPI_CONSTANTS #   
  b       OPEQ+MPI_CONSTANTS '   l  T   a   OPEQ%LHS+MPI_CONSTANTS '   �  T   a   OPEQ%RHS+MPI_CONSTANTS &     b       INFONEQ+MPI_CONSTANTS *   v  V   a   INFONEQ%LHS+MPI_CONSTANTS *   �  V   a   INFONEQ%RHS+MPI_CONSTANTS )   "  b       REQUESTNEQ+MPI_CONSTANTS -   �  Y   a   REQUESTNEQ%LHS+MPI_CONSTANTS -   �  Y   a   REQUESTNEQ%RHS+MPI_CONSTANTS &   6  b       COMMNEQ+MPI_CONSTANTS *   �  V   a   COMMNEQ%LHS+MPI_CONSTANTS *   �  V   a   COMMNEQ%RHS+MPI_CONSTANTS %   D  b       WINNEQ+MPI_CONSTANTS )   �  U   a   WINNEQ%LHS+MPI_CONSTANTS )   �  U   a   WINNEQ%RHS+MPI_CONSTANTS &   P  b       FILENEQ+MPI_CONSTANTS *   �  V   a   FILENEQ%LHS+MPI_CONSTANTS *     V   a   FILENEQ%RHS+MPI_CONSTANTS )   ^  b       MESSAGENEQ+MPI_CONSTANTS -   �  Y   a   MESSAGENEQ%LHS+MPI_CONSTANTS -     Y   a   MESSAGENEQ%RHS+MPI_CONSTANTS ,   r  b       ERRHANDLERNEQ+MPI_CONSTANTS 0   �  \   a   ERRHANDLERNEQ%LHS+MPI_CONSTANTS 0   0  \   a   ERRHANDLERNEQ%RHS+MPI_CONSTANTS '   �  b       GROUPNEQ+MPI_CONSTANTS +   �  W   a   GROUPNEQ%LHS+MPI_CONSTANTS +   E  W   a   GROUPNEQ%RHS+MPI_CONSTANTS *   �  b       DATATYPENEQ+MPI_CONSTANTS .   �  Z   a   DATATYPENEQ%LHS+MPI_CONSTANTS .   X  Z   a   DATATYPENEQ%RHS+MPI_CONSTANTS $   �  b       OPNEQ+MPI_CONSTANTS (     T   a   OPNEQ%LHS+MPI_CONSTANTS (   h  T   a   OPNEQ%RHS+MPI_CONSTANTS    �  �       DELF    H  �       ANEAR_NODE    �  �       BNEAR_NODE    `  �       YUADIST    �  �       YUBDIST    x  �       XUADIST      �       XUBDIST    �  �       SETUP_TVD #   "  @     SETUP_TVD%NCV+LIMS !   b  @     SETUP_TVD%M+LIMS %   �  ]       MPI_OP+MPI_CONSTANTS -   �  H   a   MPI_OP%MPI_VAL+MPI_CONSTANTS +   G  ]       MPI_DATATYPE+MPI_CONSTANTS 3   �  H   a   MPI_DATATYPE%MPI_VAL+MPI_CONSTANTS (   �  ]       MPI_GROUP+MPI_CONSTANTS 0   I  H   a   MPI_GROUP%MPI_VAL+MPI_CONSTANTS -   �  ]       MPI_ERRHANDLER+MPI_CONSTANTS 5   �  H   a   MPI_ERRHANDLER%MPI_VAL+MPI_CONSTANTS *   6  ]       MPI_MESSAGE+MPI_CONSTANTS 2   �  H   a   MPI_MESSAGE%MPI_VAL+MPI_CONSTANTS '   �  ]       MPI_FILE+MPI_CONSTANTS /   8  H   a   MPI_FILE%MPI_VAL+MPI_CONSTANTS &   �  ]       MPI_WIN+MPI_CONSTANTS .   �  H   a   MPI_WIN%MPI_VAL+MPI_CONSTANTS '   %   ]       MPI_COMM+MPI_CONSTANTS /   �   H   a   MPI_COMM%MPI_VAL+MPI_CONSTANTS *   �   ]       MPI_REQUEST+MPI_CONSTANTS 2   '!  H   a   MPI_REQUEST%MPI_VAL+MPI_CONSTANTS '   o!  ]       MPI_INFO+MPI_CONSTANTS /   �!  H   a   MPI_INFO%MPI_VAL+MPI_CONSTANTS 