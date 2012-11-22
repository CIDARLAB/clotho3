/*
 Copyright (c) 2009 The Regents of the University of California.
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
 ENHANCEMENTS, OR MODIFICATIONS..
 */

package org.clothocore.util.buttons;

import java.awt.AlphaComposite;
import java.awt.Graphics;
import java.awt.Graphics2D;
import javax.swing.ImageIcon;
import javax.swing.JButton;

/**
 *
 * @author J. Christopher Anderson
 */
public class TransparentButton extends JButton {


public TransparentButton(ImageIcon icon) {
    super("");
    this.setContentAreaFilled(false);
    setBorder(null);
    _myIcon = icon;
    setIcon(_myIcon);
    setTransparency(exitAlpha);
    addMouseListener(new java.awt.event.MouseAdapter() {
        @Override
        public void mouseEntered(java.awt.event.MouseEvent evt) {
            transparentButtonMouseEntered(evt);}
        @Override
        public void mouseExited(java.awt.event.MouseEvent evt) {
            transparentButtonMouseExited(evt);}
    });
}

        @Override
	public void paint(Graphics g) {
	    Graphics2D g2 = (Graphics2D) g.create();
	    g2.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, transparency));
	    super.paint(g2);
	    g2.dispose();
	}

    private void transparentButtonMouseEntered(java.awt.event.MouseEvent evt) {
        setTransparency(enterAlpha);
    }

    private void transparentButtonMouseExited(java.awt.event.MouseEvent evt) {
        setTransparency(exitAlpha);
    }

    private void setTransparency(float newalpha) {
        transparency = newalpha;
    }

    /* SETTERS
     * */

    @Override
     public void setLocation(int x, int y) {
        setBounds(x, y, _myIcon.getIconWidth(), _myIcon.getIconHeight());
     }

     public void setExitAlpha(float newalpha) {
         exitAlpha = newalpha;
         setTransparency(exitAlpha);
     }

     public void reset() {
         setTransparency(exitAlpha);
         repaint();
     }

     public void setEnterAlpha(float newalpha) {
         enterAlpha = newalpha;
     }

     public void changeIcon(ImageIcon newIcon) {
        _myIcon = newIcon;
        setIcon(_myIcon);
     }

    /* PUTTERS
     * */

    /* GETTERS
     * */

/*-----------------
     variables
 -----------------*/

    private float transparency = 0.5f;
    private float exitAlpha = 0.2f;
    private float enterAlpha = 1.0f;
    private ImageIcon _myIcon;
}
