/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.clothocore.util.buttons;

import java.awt.AlphaComposite;
import java.awt.Graphics;
import java.awt.Graphics2D;
import javax.swing.JButton;

/**
 *
 * @author jcanderson
 */
public class GlassButton  extends JButton {
	public GlassButton(String text) {
	    super(text);
	    setOpaque(false);
	}

	public GlassButton(float trans) {
	    super("");
            transparency = trans;
	    setOpaque(false);
	}

    @Override
    public void paint(Graphics g) {
        Graphics2D g2 = (Graphics2D) g.create();
        g2.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, transparency));
        super.paint(g2);
        g2.dispose();
    }

    private float transparency = 0.0f;

}
