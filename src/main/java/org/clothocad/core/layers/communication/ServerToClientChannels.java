/*
Copyright (c) 2010 The Regents of the University of California.
All rights reserved.
Permission is hereby granted, without written agreement and without
license or royalty fees, to use, copy, modify, and distribute this
software and its documentation for any purpose, provided that the above
copyright notice and the following two paragraphs appear in all copies
of this software.

IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.

THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE
PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
ENHANCEMENTS, OR MODIFICATIONS.
 */
package org.clothocad.core.layers.communication;

/**
 * Needs to match up with `web/scripts/clotholib/librecv.js`
 *
 * @author Kelvin Li
 */
public enum ServerToClientChannels {
    //Command bar functionality
    INTERNAL_CHANNEL,
    showCommandResult,     //For displaying the potential command results of command bar
    showQueryCompletions,  //The autocomplete, one-character-at-a-time or in response to disambiguation
    showTabList,           //Send links, maybe a screenshot, to display the workspaces a user has configured
    login,                 //JCA:  SOMETHING ABOUT LOGIN COMMUNICATIONS
    logout,                //JCA:  SOMETHING ABOUT THE LOGOUT COMMUNICATIONS
    
    commandList            //Specific HTML-related commands for manipulating the view
    
//    showWidget,
//    updateWidget,
//    collect
//    addPage,
//    clearAuthKey,
//    removePage,
//    removeWidget,
//    setAuthKey,
};
